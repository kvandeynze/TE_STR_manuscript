import tdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import subprocess
import os
import multiprocessing
import re

# Load the data and output the allele table
def load_data(database_path, exclude_locus_ids):
    # Connect to the database
    db = tdb.load_tdb(database_path, 
                    lfilters=[( 'LocusID', 'not in', exclude_locus_ids)])
    
    return db['allele']

#function to write the sequence on each allele to a fasta file with name corresponding to LocusID:allele_number
def write_fasta(data, output_path):
    with open(output_path, 'w') as f:
        for i, row in data.iterrows():
            f.write(f'>{row.LocusID}:{row.allele_number}\n')
            f.write(f'{row.sequence}\n')
    return

#function to run TRF on the fasta file using subprocess
def run_trf(fasta_path):
    trf_path = '/data/software/TRF/4.09.1/trf409.linux64'  # Full path to the trf executable
    trf_command = [trf_path, fasta_path, '2', '5', '5', '80', '10', '10', '500', '-d', '-h']
    # trf_command = ['trf', fasta_path, '2', '5', '5', '80', '10', '10', '500', '-d', '-h']
    result = subprocess.run(trf_command, capture_output=True, text=True)
    
    # Construct the .dat file name
    fasta_filename = os.path.basename(fasta_path)
    dat_file = fasta_filename + '.2.5.5.80.10.10.500.dat'
    # Check if the .dat file exists
    if not os.path.exists(dat_file):
        raise FileNotFoundError(f"TRF did not produce the expected .dat file: {dat_file}")
    return dat_file, fasta_path
def run_trf_multiprocessing(fasta_files):
    with multiprocessing.Pool() as pool:
        dat_files_and_fasta_files = pool.map(run_trf, fasta_files)
    return dat_files_and_fasta_files

#function to parse the TRF output and return the percent of the input sequence that is a tandem repeat -- curently only checks % detected as a TR not the specific motifs -- we can see how this works initially
def parse_trf_output(dat_file, fasta_file):
    # Read the input sequences from the fasta file
    sequences = {}
    with open(fasta_file, 'r') as f:
        sequence_name = None
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name:
                    sequences[sequence_name] = ''.join(sequence)
                sequence_name = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if sequence_name:
            sequences[sequence_name] = ''.join(sequence)
    
    # Initialize a dictionary to store the percentages
    percentages = {name: 0 for name in sequences.keys()}
    
    # Regular expression to match lines starting with a numeric value
    numeric_re = re.compile(r'^\d+')
    
    # Parse the TRF .dat output file
    with open(dat_file, 'r') as f:
        current_sequence = None
        intervals = []
        for line in f:
            line = line.strip()
            if line.startswith('Sequence:'):
                if current_sequence and intervals:
                    # Merge overlapping intervals
                    merged_intervals = []
                    for start, end in sorted(intervals):
                        if merged_intervals and merged_intervals[-1][1] >= start - 1:
                            merged_intervals[-1] = (merged_intervals[-1][0], max(merged_intervals[-1][1], end))
                        else:
                            merged_intervals.append((start, end))
                    
                    # Calculate the total length of tandem repeats
                    total_tr_length = sum(end - start + 1 for start, end in merged_intervals)
                    
                    # Calculate the percentage of the sequence that is tandem repeats
                    sequence_length = len(sequences[current_sequence])
                    percentages[current_sequence] = (total_tr_length / sequence_length) * 100
                
                current_sequence = line.split(' ')[1].strip()
                intervals = []
            elif numeric_re.match(line):
                trf_line = line.split()
                start = int(trf_line[0])
                end = int(trf_line[1])
                intervals.append((start, end))
        
        # Handle the last sequence
        if current_sequence and intervals:
            # Merge overlapping intervals
            merged_intervals = []
            for start, end in sorted(intervals):
                if merged_intervals and merged_intervals[-1][1] >= start - 1:
                    merged_intervals[-1] = (merged_intervals[-1][0], max(merged_intervals[-1][1], end))
                else:
                    merged_intervals.append((start, end))
            
            # Calculate the total length of tandem repeats
            total_tr_length = sum(end - start + 1 for start, end in merged_intervals)
            
            # Calculate the percentage of the sequence that is tandem repeats
            sequence_length = len(sequences[current_sequence])
            percentages[current_sequence] = (total_tr_length / sequence_length) * 100
    
    return percentages

def merge_intervals(intervals):
    merged_intervals = []
    for start, end in sorted(intervals):
        if merged_intervals and merged_intervals[-1][1] >= start - 1:
            merged_intervals[-1] = (merged_intervals[-1][0], max(merged_intervals[-1][1], end))
        else:
            merged_intervals.append((start, end))
    return merged_intervals

def parse_trf_output_periods(dat_file, fasta_file, periods_locusid):
        #TODO make sure this works, add an additional period column that is the percent of the consensnsus size that is the period that matches the input period since sometimes TRF likes to get slightly larger consensus from the longest or shortest possible period
        # Read the input sequences from the fasta file
        sequences = {}
        with open(fasta_file, 'r') as f:
            sequence_name = None
            sequence = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if sequence_name:
                        sequences[sequence_name] = ''.join(sequence)
                    sequence_name = line[1:]
                    sequence = []
                else:
                    sequence.append(line)
            if sequence_name:
                sequences[sequence_name] = ''.join(sequence)
        
        # Initialize a dictionary to store the periods
        periods = {name: {'longest': 0, 'smallest': float('inf'), 'largest': 0,"percent_match_input_period":float('-inf')} for name in sequences.keys()}
        
        # Regular expression to match lines starting with a numeric value
        numeric_re = re.compile(r'^\d+')
        
        # Parse the TRF .dat output file
        with open(dat_file, 'r') as f:
            current_sequence = None
            intervals = []
            for line in f:
                line = line.strip()
                if line.startswith('Sequence:'):
                    # given all intervals where the period matches the input period, merge them and calculate the percentage of the sequence that is the input period
                    if intervals and current_sequence:
                        merged_intervals = merge_intervals(intervals)
                        total_tr_length = sum(end - start + 1 for start, end in merged_intervals)
                        sequence_length = len(sequences[current_sequence])
                        periods[current_sequence]["percent_match_input_period"] = (total_tr_length / sequence_length) * 100
                    current_sequence = line.split(' ')[1].strip()
                    intervals = []
                    try:
                        curr_input_period = periods_locusid[float(current_sequence.split(":")[0])] # get period from the locusid dictionary
                        # print(f"current input period: {curr_input_period}")
                    except:
                        # print(f"Warning: LocusID {current_sequence.split(':')[0]} not found in periods_locusid file")
                        curr_input_period = 0
                        curr_percent_match_input_period = 0 # keep this to make this cumlative
                elif numeric_re.match(line):
                    #initialize all period values to the first period value because if there is only 1 entry it is the longest, smallest, and largest
                    if periods[current_sequence]['longest'] == 0:
                        trf_line = line.split()
                        period = int(trf_line[2])
                        periods[current_sequence]['longest'] = period
                        periods[current_sequence]['smallest'] = period
                        periods[current_sequence]['largest'] = period
                        curr_consensus_size = len(trf_line[-1])

                    trf_line = line.split()
                    period = int(trf_line[2])
                    consensus_size = len(trf_line[-1])
                    start = int(trf_line[0])
                    end = int(trf_line[1])

                    if period == curr_input_period:
                        intervals.append((start, end)) # intervals is just intervals for matching period
                        # FIXME this curently is incorrect because we are not taking into consideration overlappoing TRF results of the same period
                        # record the percentage of the sequence that can be identified as the input period
                        # periods[current_sequence]["percent_match_input_period"] = ((consensus_size + curr_percent_match_input_period) / len(sequences[current_sequence])) * 100 # the only scenario I see this not working well is if there are multiple stretches of the same period in the sequence
                        # curr_percent_match_input_period = consensus_size + curr_percent_match_input_period # update the current size of the consensus to keep track of the total size of the consensus matching the input period
                    if period > periods[current_sequence]['longest']:
                        periods[current_sequence]['longest'] = period
                    if period < periods[current_sequence]['smallest']:
                        periods[current_sequence]['smallest'] = period
                    if consensus_size >= curr_consensus_size: # we want to keep the shortest period per consensus size
                        if consensus_size == curr_consensus_size:
                            if period < periods[current_sequence]['largest']: # leave the period be if it is the same as the current largest and the current period is greater than the current largest
                                periods[current_sequence]['largest'] = period
                        else:
                            periods[current_sequence]['largest'] = period
                        curr_consensus_size = consensus_size # updated so we aren't always overriding the last largest by comparing to period, need to update to the consensus size

            # Replace 'inf' with 0 if no repeats were found
            for key in periods:
                if periods[key]['smallest'] == float('inf'):
                    periods[key]['smallest'] = periods[key]['longest'] # there was a bug here that was setting the last current sequence longest for the multiprocess child process to all smallest values when not available

            # deal with last sequence for merging intervals
            if intervals and current_sequence:
                merged_intervals = merge_intervals(intervals)
                total_tr_length = sum(end - start + 1 for start, end in merged_intervals)
                sequence_length = len(sequences[current_sequence])
                periods[current_sequence]["percent_match_input_period"] = (total_tr_length / sequence_length) * 100
        
        return periods
# Function to parse TRF outputs in parallel and combine the results
def parse_trf_output_multiprocessing(dat_files_and_fasta_files):
    with multiprocessing.Pool() as pool:
        results = pool.starmap(parse_trf_output, dat_files_and_fasta_files)
    
    # Combine the results into a single dictionary
    combined_percentages = {}
    for result in results:
        combined_percentages.update(result)
    
    return combined_percentages
def parse_trf_output_periods_multiprocessing(dat_files_and_fasta_files):
    with multiprocessing.Pool() as pool:
        results = pool.starmap(parse_trf_output_periods, dat_files_and_fasta_files)
    
    # Combine the results into a single dictionary
    periods = {}
    for result in results:
        periods.update(result)
    
    return periods
def read_integers_from_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read().strip()
        
        # Remove the square brackets
        content = content.strip('[]')
        
        # Split the string by commas and convert to integers
        integer_list = [int(num.strip()) for num in content.split(',')]
        
    return integer_list
def read_fasta(fileObject):
  '''
  Generator function to read in reads from fasta
  '''
  header = ''
  seq = ''
  # skip any useless leading information
  for line in fileObject:
    if line.startswith('>'):
      header = line.strip()
      break
  for line in fileObject:
    if line.startswith('>'):
      if 'N' in seq.upper():
        yield header, seq, False
      else:
        yield header, seq, True
      header = line.strip()
      seq = ''
    else:
      seq += line.strip()
  if header:
    if 'N' in seq.upper():
      yield header, seq, False
    else:
      yield header, seq, True
def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('database_path', type=str, help='Path to the tdb directory')
    parser.add_argument('exclude_locus_ids', type=str, help='File path to list of LocusIDs to exclude')
    parser.add_argument('output_path', type=str, help='Path to save the output')

    args = parser.parse_args()

    # read the list from the file 
    exclude_locus_ids = read_integers_from_file(args.exclude_locus_ids)
    # data = load_data(args.database_path, exclude_locus_ids)

    # write the alleles to a fasta file to run trf on consensus sequences
    # write_fasta(data, args.output_path +'_alleles.fasta') # for now commenting this out since I already ran this, TODO update this to split the fasta files into smaller chunks and run TRF in parallel

    # Run TRF on the fasta file
    # trf_output_file = run_trf(args.output_path + '_alleles.fasta')

    # Parse the TRF output
    # percentages = parse_trf_output(trf_output_file, args.output_path + '_alleles.fasta')

    #multiprocessing version
    fasta_files = [os.path.join('/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC', f) 
                   for f in os.listdir('/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC') 
                   if f.startswith('all_HPRC_100_with_AP_allelesa') and f.endswith('.fa')]
    # dat_files_and_fasta_files = run_trf_multiprocessing(fasta_files)
    print("dat files produced! moving on to parsing output...")
    # Reconstruct the dat_files_and_fasta_files by pairing fasta files with corresponding dat files
    # read periods per locusid
    periods_locusid = pd.read_csv('/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/merged_starting_MEIs/MEI_and_non_MEI_STRs.period.locusid.txt', sep=' ')
    periods_locusid.columns = ['period','LocusID']
    periods_locusid["period"] = periods_locusid["period"].astype(float)
    periods_locusid["LocusID"] = periods_locusid["LocusID"].astype(float)

    # periods locusid is a dataframe that contains the period for each locusid. each locusid is the first part of the sequence name in the fasta file ie 148:1 the locus id would be 148
    # we want to suvset the periods_locusid dataframe to only include the locusids that are in the fasta files per fasta file dat pair
    dat_files_and_fasta_files = []
    for fasta_file in fasta_files:
        fasta_filename = os.path.basename(fasta_file)
        dat_file = fasta_filename + '.2.5.5.80.10.10.500.dat'
        # get the locusid from the fasta file and subset the periods_locusid dataframe to only include the locusids
        # get all locusids from fasta file headers
        locusids = set([float(header.split('>')[1].split(":")[0]) for header,seq,indicator in read_fasta(open(fasta_file))])
        # subset the periods_locusid dataframe to only include the locusids in the fasta file
        periods_locusid_subset = periods_locusid[periods_locusid['LocusID'].isin(locusids)]
        periods_dict = periods_locusid_subset.set_index('LocusID')['period'].to_dict()

        if os.path.exists(dat_file):
            dat_files_and_fasta_files.append((dat_file, fasta_file, periods_dict))
        else:
            print(f"Warning: Expected .dat file not found for {fasta_file}")

    print("dat files and fasta files paired! moving on to parsing output...")
    # percentages = parse_trf_output_multiprocessing(dat_files_and_fasta_files) # this is actually REALLY fast

    periods = parse_trf_output_periods_multiprocessing(dat_files_and_fasta_files) # this is actually REALLY fast
    
    # Save the data to the output path
    # percentages_df = pd.DataFrame(list(percentages.items()), columns=['Allele', 'Percentage'])
    # percentages_df.to_csv(args.output_path + "_percent_TR.tsv", index=False, sep='\t')
    periods_df = pd.DataFrame.from_dict(periods, orient='index')
    periods_df.to_csv(args.output_path + "_periods_revised_input_match_2.tsv", sep='\t')

if __name__ == '__main__':
    main()