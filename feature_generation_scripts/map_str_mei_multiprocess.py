import matplotlib.pyplot as plt
from Bio import Phylo
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
import pandas as pd
import mappy as mp
import requests
import multiprocessing

from subset_alu_seqs import *
from analysis_utils import *
import argparse
import pickle
from functools import partial
import os

def fetch_consensus_sequences(family_name):
    """
    Fetches the consensus sequence for a given family from DFAM.
    
    Parameters:
        family_name (str): The name of the family to query.
    
    Returns:
        str: The consensus sequence for the family.
    """
    url = "https://dfam.org/api/families"

    params = {
        "name": family_name,
        "clade": "Primates",
        "clade_relatives": "both"

    }

    response = requests.get(url, params=params)
    results = response.json()["results"]

    all_consensus = {}
    for element in results:
        # get the consensus from the id request
        url = f"https://dfam.org/api/families/{element['accession']}"
        response = requests.get(url)
        all_consensus[element["name"]] = response.json()["consensus_sequence"]

    # all consensus sequences have 30bp polyA tail, trim them all
    for key in all_consensus:
        all_consensus[key] = all_consensus[key][:-30]

    return all_consensus

def count_gaps_up_to_position(sequence, position):
    """
    This is an updated function that is aware that the position is going to be shifted due to the
    gaps we are counting. That is, this function will count the gaps up until position not including the gaps.
    """
    gap_count = 0
    curr_pos = 0
    while curr_pos < position:
        if sequence[curr_pos] == "-":
            gap_count += 1
            position += 1  # Shift the position forward since we are counting a gap
        curr_pos += 1
    return gap_count

def adjust_positions_for_multisequence(family_name, positions,sequences):
    """
    Adjusts the positions for a given family in the alignment by accounting for gaps.
    This is the inverse of adjusting the positions from the overall consensus alignment to the family alignment.
    Here we are taking positions from the family alignment and adjusting them to the overall consensus alignment
    so we want to add the gaps instead of subtracting them.
    
    Parameters:
        alignment_file (str): Path to the ClustalW alignment file.
        family_name (str): The name of the family to query.
        positions (list): List of 1-based positions to adjust.
    
    Returns:
        list: Adjusted positions accounting for gaps.
    """
    
    # Get the sequence for the specified family
    if family_name not in sequences:
        raise ValueError(f"Family '{family_name}' not found in the alignment.")
    family_sequence = sequences[family_name]
    
    # Adjust positions
    adjusted_positions = []
    for pos in positions:
        # Convert 1-based position to 0-based index
        zero_based_pos = pos - 1
        # Count gaps up to the position

        if zero_based_pos < 0:
            print(f"Zero-based position {zero_based_pos} is negative, counting no gaps.")
            gap_count = 0
        elif zero_based_pos >= len(family_sequence.replace("-",""))-1:
            print(f"Zero-based position {zero_based_pos} is beyond the sequence length, counting all gaps.")
            # gap_count = count_gaps_up_to_position(family_sequence, len(family_sequence) - 1)
            # print(f"Gap count: {gap_count}")
            gap_count = family_sequence.count("-")
        else:
            gap_count = count_gaps_up_to_position(family_sequence, zero_based_pos)
        # Adjust the position
        adjusted_positions.append(zero_based_pos + gap_count)
    return adjusted_positions
def adjust_relative_positions(query_aln, str_relative_start):
    # FIXME this function counts the gap UP UNTIL the relative start position
    # For STRs that are upstream of the MEI, the relative start position is negative so we are always adding all gaps to this which is incorrect
    # I beleive we can just return the original relative start position for these cases but I need to think about if that is correct or not
    if str_relative_start < 0:
        # check how many gaps are in the query alignment before the start position according to
        # the family alignment (leading gaps)

        # FIXME if relative start is negative, we need to adjust based on the resverse complement
        # print(f"Relative start position {str_relative_start} is negative, returning gaps before the start position.")
        gap_start = 0
        for i in range(len(query_aln)):
            if query_aln[i] == "-":
                gap_start += 1
            else:
            # Stop counting once we hit the first non-gap character
                break
        new_relative_start = str_relative_start - gap_start

        return new_relative_start
    # calculate the new relative start and end positions by accounting for the gaps up until the relative start
    trailing_query_gaps = get_trailing_gaps(query_aln)
    gap_start = 0
    for i in range(len(query_aln)-trailing_query_gaps): #FIXME changed this to only iterate up until the trailing
        if query_aln[i] == "-":
            gap_start += 1
        # keep track of how many non-gap characters we have seen to see if we have reached the relative start
        if i - gap_start == str_relative_start:
            break
    new_relative_start = str_relative_start + gap_start

    return new_relative_start
# Count the number of "-" characters at the start of the reversed alignment before a non-gap character
def get_trailing_gaps(family_align):
    count = 0
    for c in family_align[::-1]:
        if c == "-":
            count += 1
        else:
            break
    return count

def get_consensus_position(family_align, query_aln, str_relative_start):
    new_relative_start = adjust_relative_positions(query_aln, str_relative_start)
    # print(f"New relative start: {new_relative_start}")

    # FIXME This function assumes that the relative start is positive and thus does not account for upstream positions 
    # here assume that we don'tneed to account for gaps if the relative start is negative, THIS MAY BE NAIIVE AND INCORRECT
    # if new_relative_start < 0:
    #     print(f"Relative start position {new_relative_start} is negative, returning original position.")
    #     return new_relative_start

    # get the position with respect to the family alignment
    family_position = 0
    # only go up until the end of the family consensus alignment based on trailing gaps, ie don't iterate into the gaps at the end of the alignment
    # get the number of trailing gaps in the family alignment
    
    # NOTE this will only ever return for internal start positions, not for upstream or poly(A)
    trailing_gaps = get_trailing_gaps(family_align)
    # print(f"Trailing gaps: {trailing_gaps}")
    for i in range(len(family_align) - trailing_gaps):
        if family_position == new_relative_start:
            # get the number of gaps up to this point and subtract them from the index
            final = i - family_align[:i].count("-")
            return final
        family_position += 1
    

    # FIXME I need to refactor this for clarity, this is a bit convoluted -- final should just be new_relative_start minus the number of gaps up to the poly(A)
    # we still need the number of trailing gaps to acount for only counting internal gaps since the consensus doesn't include the Poly(A) tail
    # however, the new_relative_start is the position relative to both alignments so its relative start is just its position in the alignment minus the non-trailing gaps
    # since the new_relative_position is the position relative to the end of the consensus including the poly(A) sequence included (or not included) in the reference
    # Alu sequence

    # NOTE added after debugging for clarity
    consensus_end = len(family_align) - trailing_gaps
    gaps_up_to_polyA = family_align[:consensus_end].count("-")
    
    # account for upstream positions -- the sign reflects which direction the relative start is in with
    # respect to the start of the family alignment, the magnitude reflects the distance from the start
    # thus we need to account for if the gaps are making the relative start more negative or more positive
    if new_relative_start < 0:
        # we don't care about any gaps after the start position since it is upstream of the MEI
        # get the number of leading gaps in the family alignment
        # print(f"Relative start position {new_relative_start} is negative, adjusting for gaps before the start position.")
        gap_start = 0
        for i in range(len(family_align)):
            if family_align[i] == "-":
                gap_start += 1
            else:
            # Stop counting once we hit the first non-gap character
                break
        final = new_relative_start + gap_start
    else:
        final = new_relative_start - gaps_up_to_polyA

    # print(f"Final position: {final}")

    return final

# function to convert reverse oriented relative start and end positions to the forward oriented version of the sequence 
def convert_reverse_relative_positions_to_forward(str_relative_start, str_relative_end, seq_length):
    # the old relative end is now the new relative start, and the old relative start is now the new relative end
    # FIXME seq_length - 1 - str_relative_end, seq_length - 1 - str_relative_start (accounts for length)
    return seq_length - str_relative_end - 1 , seq_length - str_relative_start - 1

# function that takes the mei_id, Alu_subfamily, orientation, start, end, start_mei, and end_mei columns as well as the consensus dictionary to pull from for a given row and returns the start of the STR relative to the MEI consensus sequence
# FIXME be on the look out for off by one errors in the start and end positions of the strs
def get_consensus_relative_start(row, all_alu_consensus,family_alignments, seq):
    # get the Alu subfamily from the row
    alu_subfamily = row["mei_subfamily"]
    # get the orientation from the row
    orientation = row["orientation"]
    # get the start and end positions of the STR in the Alu sequence
    str_start = row["start"]
    str_end = row["end"]
    # get the start and end positions of the MEI in the Alu sequence
    start_mei = row["start_mei"]
    end_mei = row["end_mei"]

    # get the Alu consensus sequence for this subfamily
    try:
        family_seq = all_alu_consensus[alu_subfamily]
    except KeyError:
        print(f"Alu consensus sequence for {alu_subfamily} not found, returning None")
        return None

    # get the relative start and end positions of the STR in the Alu sequence
    # FIXME make sure the inputs have accounted for the off by 1
    str_relative_start = str_start - start_mei
    str_relative_end = str_end - start_mei - 1 # minus one to account for exclusive end
    # print(f"Relative start: {str_relative_start}, Relative end: {str_relative_end}")
    if orientation == "-":
        seq_len = end_mei - start_mei #FIXME end_mei - start_mei (account for the buffer -1 in input start)
        str_relative_start, str_relative_end = convert_reverse_relative_positions_to_forward(str_relative_start, str_relative_end, seq_len)
        # seq = mp.revcomp(seq) # get the reverse complement of the sequence
        # print(f"Converted relative start: {str_relative_start}, relative end: {str_relative_end}")

        # FIXME assume that if the relative start is negative then we don't need to adjust anything since it is upstream of the MEI and the consensus sequence is not affected by changes we make to the sequence
    # if str_relative_start < 0:
    #     print(f"Relative start position {str_relative_start} is negative, returning original position.")
    #     return str_relative_start

        #get the Alu sequence from the aligner using the mei_id -- take this as an input instead
    # seq = aligner.seq(row["mei_id"])
    # print(f"The extracted forward oriented STR sequence is: {seq[str_relative_start:str_relative_end]}")
    # print(f'motif was: {row["str_motif"]}')

    # if seq is empty, the Alu sequence is not found in our full length Alu set, return None
    if not seq:
        # print(f"Alu sequence for {row['mei_id']} not found in full length Alu set, returning None")
        return None

    # align the family sequence to the Alu sequence
    family_align, query_aln = align_family_sequences(family_seq, seq)
    # print(family_align)
    # print(query_aln)

    curr_consensus_start = get_consensus_position(family_align, query_aln, str_relative_start)
    # adjust for the gaps in the family alignments to cast all positions to the multisequence alignment of the families -- ASSUMES 1 BASED INDEXES FOR THE POSITIONS
    final_position = adjust_positions_for_multisequence(alu_subfamily, [curr_consensus_start + 1], family_alignments) # FIXME this converts positions TO the family sequence not back, we want to ADD the gaps
    # print(f"Final position in the multisequence alignment: {final_position[0]}")
    return final_position[0]

def initializer(global_all_consensus, global_family_alignments):
    global all_consensus, family_alignments
    all_consensus = global_all_consensus
    family_alignments = global_family_alignments
def process_row_with_seq(args):
    row, seq, all_consensus, family_alignments = args
    return get_consensus_relative_start(row, all_consensus, family_alignments, seq)

def main():

    parser = argparse.ArgumentParser(description="Map STRs to MEI consensus positions.")
    parser.add_argument("--input", required=True, help="Input tab delimited file with STR and MEI start and end positions, NOTE: assumes no header is included")
    parser.add_argument("--family", required=True, help="MEI family name to query consensus sequences from DFAM")
    parser.add_argument("--family_msa", required=True, help="Alignment file with family alignments in ClustalW format")
    parser.add_argument("--fasta", required=True, help="FASTA file with MEI sequences from the reference genome")
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPUs to use for multiprocessing (default: 1)")
    args = parser.parse_args()

    # Load input data -- assume no header
    position_data = pd.read_csv(args.input, sep="\t", header=None)
    try:
        position_data.columns = ["chrom", "start", "end", "LocusID", "chrom_mei", "start_mei","end_mei","mei_id","mei_subfamily","orientation","dist"]
    except ValueError:
        print("Input file does not have the expected number of columns. Please check the input format. Expects output from bedtools closest -d")
        print("Check get_1kb_alu_strs.sh for the expected input format.")
        print("Expected columns: chrom, start, end, LocusID, chrom_mei, start_mei, end_mei, mei_id, mei_subfamily, orientation")
        exit(1)
    print(f"Input data loaded with {len(position_data)} rows.")

    # load the family consensus sequences from DFAM
    all_consensus = fetch_consensus_sequences(args.family)

    # load the family alignments from the ClustalW alignment file
    family_alignments = parse_clustalw_alignment(args.family_msa)

    # load the full length Alu sequences from the FASTA file to the aligner
    aligner = mp.Aligner(args.fasta) # TODO make sure passing the aligner to a multiprocessing pool works correctly
    if not aligner:
        print("Error: could not load the FASTA file with full length Alu sequences.")
        exit(1)
    print(f"Aligner loaded.")

    # Prepare arguments for multiprocessing
    # instead of passing the aligner, get the sequences from the aligner and pass that with the row
    # This is to avoid issues with pickling the aligner object
    # Prepare a list of (row, sequence, all_consensus, family_alignments) tuples to process
    row_seq_tuples = [
        (row, aligner.seq(row["mei_id"]), all_consensus, family_alignments)
        for _, row in position_data.iterrows()
    ]

    # Update process_row to accept (row, seq, all_consensus, family_alignments) tuple


    with multiprocessing.Pool(processes=args.cpus) as pool:
        results = pool.map(process_row_with_seq, row_seq_tuples)
    # Collect results
    position_data["consensus_relative_start"] = results
    # Save the results to a new file
    base, ext = os.path.splitext(args.input)
    output_file = base + "_consensus_relative_start.txt"
    position_data.to_csv(output_file, sep="\t", index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
