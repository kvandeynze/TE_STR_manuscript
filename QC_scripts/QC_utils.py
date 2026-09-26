import tdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sb
from scipy.stats import pearsonr
import matplotlib.colors as mcolors
import functools
# import from the directory above this one
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd().parent))
from analysis_utils import TDBUtils
import logging

class tdb_qc:
    # TODO consider changing this to not being static, save data as an attribute of the class instead

    @staticmethod
    def log_changes(func):
        @functools.wraps(func)
        def wrapper_log_changes(*args, **kwargs):
            df = args[0] 
            
            # check if LocusID:allele_number is in columns
            if 'LocusID:allele_number' not in df.columns and "LocusID" in df.columns and "allele_number" in df.columns:
                df['LocusID:allele_number'] = df['LocusID'].astype(str) + ":" + df['allele_number'].astype(str)

            initial_alleles = df['LocusID:allele_number'].nunique() if 'LocusID:allele_number' in df.columns else 0
            initial_loci = df['LocusID'].nunique() if 'LocusID' in df.columns else 0
            
            result = func(*args, **kwargs)

            if isinstance(result, tuple):
                result_df = result[0]
            else:
                result_df = result
            
            # check if LocusID:allele_number is in columns
            if 'LocusID:allele_number' not in result_df.columns and "LocusID" in result_df.columns and "allele_number" in result_df.columns:
                result_df['LocusID:allele_number'] = result_df['LocusID'].astype(str) + ":" + result_df['allele_number'].astype(str)

            final_alleles = result_df['LocusID:allele_number'].nunique() if 'allele_number' in result_df.columns else 0
            final_loci = result_df['LocusID'].nunique() if 'LocusID' in result_df.columns else 0
            
            print(f"Function {func.__name__}:")
            print(f"  Alleles: {initial_alleles} -> {final_alleles}")
            print(f"  Loci: {initial_loci} -> {final_loci}")
            
            return result
        return wrapper_log_changes
    
    @staticmethod
    def load_period_data(period_loci_path, period_stats_path):
        '''
        Load the period data and return a pandas dataframe 
        '''
        period_locusids = pd.read_csv(period_loci_path, sep=' ', header=None)
        period_locusids.columns = ["period", "LocusID"]
        period_stats = pd.read_csv(period_stats_path, sep="\t")
        period_stats.columns = ["LocusID:allele", "longest_period","shortest_period","most_abundant_period","percent_input_match"]
        #split the LocusID:allele column into LocusID and allele columns
        period_stats[["LocusID","allele_number"]] = period_stats["LocusID:allele"].str.split(":", expand=True)
        period_stats["LocusID"] = period_stats['LocusID'].astype(int)
        period_stats["allele_number"] = period_stats['allele_number'].astype(int)

        # merge period_stats with period_locusids
        period_stats = pd.merge(period_stats, period_locusids, on="LocusID", how = "left")

        # print("There are ", period_stats.shape[0], " alleles in the period_stats dataframe")
        # print("There are ", period_stats["LocusID"].nunique(), " loci in the period_stats dataframe")

        return period_stats
    
    @staticmethod
    def load_repeatmasker_data(repeatmasker_path):
        '''
        Load the repeatmasker data and return a pandas dataframe
        '''
        # repeatmasker inputs
        repeatmasker = pd.read_csv("/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/repeatmasker_filter/all_HPRC_100_with_AP_alleles_all.out", sep="\t", header=None)
        # drop the first row
        repeatmasker = repeatmasker.dropna()
        repeatmasker.columns = ["LocusID:allele","query_start","query_end","subfamily","family"]
        # split the LocusID:allele column
        repeatmasker["LocusID"] = repeatmasker["LocusID:allele"].str.split(":", expand=True)[0].astype(int)
        repeatmasker["allele_number"] = repeatmasker["LocusID:allele"].str.split(":", expand=True)[1].astype(int)

        # calulate the length of the element detected by repeatmasker calculate_repeatmasker_span FIXME I dont think this is necessary
        repeatmasker["element_length"] = repeatmasker["query_end"].astype(int) - repeatmasker["query_start"].astype(int)


        return repeatmasker

    @staticmethod
    @log_changes
    def combine_input_data(data, period_data, repeatmasker_data): # TODO update to only include the repeatmasker data, may not need at all
        '''
        Combine the input data with the period and repeatmasker data
        parameters:
            data: data['allele'] table, pandas dataframe with columns 'LocusID','allele_number', 'percent_input_match','percent_element', 'allele_length', and 'period'
            period_data: pandas dataframe with columns 'LocusID','allele_number', 'longest_period','shortest_period','most_abundant_period','percent_input_match'
            repeatmasker_data: pandas dataframe with columns 'LocusID','allele_number','query_start','query_end','subfamily','family','element_length'
        '''
        # merge the period data with the input data
        data = pd.merge(data, period_data, on=["LocusID","allele_number"], how = "left")
        data["LocusID:allele"] = data["LocusID"].astype(str) + ":" + data["allele_number"].astype(str)

        # merge the repeatmasker data with the input data
        data = pd.merge(data, repeatmasker_data, on=["LocusID","allele_number"], how = "left")
        data["LocusID:allele"] = data["LocusID"].astype(str) + ":" + data["allele_number"].astype(str)

        return data

    @staticmethod
    @log_changes
    def allele_QC(data):
        '''
        Filter out alleles based on their sequence content including non-TR region inclusion
        parameters:
            data: pandas dataframe with columns 'LocusID','allele_number', 'percent_input_match','percent_element', 'allele_length', and 'period'
        returns:
            cleaned_data: pandas dataframe with the same columns as data 
            to_drop: pandas dataframe with the alleles that were filtered out so we can subset when loading the data
        '''
        # filter 1: remove alleles with high percent repeatmasker annotation and low percent input match (TE inclusion)
        data_2 = data[np.where((data["allele_percent_element"] <=0) | (data["percent_input_match"] > 60),True, False)]
        # filter 2: remove alleles that are outliers in terms of percent input period (neighboring TR and flanking sequence inclusion)

        #fill -inf values with 0
        data_2['percent_input_match'] = data_2['percent_input_match'].replace(-np.inf, 0)
        # calculate the outliers using IQR
        Q1 = data_2['percent_input_match'].quantile(0.25)
        Q3 = data_2['percent_input_match'].quantile(0.75)
        IQR = Q3 - Q1

        outliers_percent_input_match = data_2[(data_2['percent_input_match'] < (Q1 - 1.5 * IQR))] # only filter out the lower outliers -- should allow 100% input match alleles to be included

        # plot the distribution of percent input match and a red line at the cutoff
        plt.hist(data_2['percent_input_match'], bins=100)
        plt.axvline(x=(Q1 - 1.5 * IQR), color='r', label='IQR cutoff')
        # add labels and title
        plt.xlabel('Percent input match')
        plt.ylabel('Frequency')
        plt.title('Distribution of percent input match')
        plt.legend()
        plt.show()

        # FIXME the final outliers should be based on ALL filtered alleles
        final_outliers = outliers_percent_input_match[~((outliers_percent_input_match.allele_length < 40)&(outliers_percent_input_match.allele_length > outliers_percent_input_match.period))] # filter these from the outlier list
        clean_data = data_2[~data_2.index.isin(final_outliers.index)] # filter these from the non-outlier list
        
        to_drop = data[~data["LocusID:allele"].isin(clean_data["LocusID:allele"])] # get the alleles that were filtered out
        
        return clean_data, to_drop

    @staticmethod
    @log_changes
    def update_AC(allele_counts, to_drop): # still good, will update what is in to_drop
        """
        Update allele counts after filtering by droppping alleles in to_drop

        Parameters:
        allele_counts: pandas DataFrame with columns 'LocusID', 'allele_number', 'AC', 'AF'
        to_drop: pandas DataFrame with columns 'LocusID', 'allele_number' returned by allele_QC

        Returns:
        clean_AC: cleaned allele_counts DataFrame
        """
        clean_AC = allele_counts[~allele_counts.set_index(['LocusID', 'allele_number']).index.isin(to_drop.set_index(['LocusID', 'allele_number']).index)]
        return clean_AC

    @staticmethod
    @log_changes
    def filter_low_support_loci(cleaned_AC, min_alleles):
        '''
        Filter out loci with less than min_alleles alleles in clean_AC and return the loci that are filtered out
        parameters:
            cleaned_AC: pandas dataframe with columns 'LocusID','allele_number', 'AC', and 'AF'
            min_alleles: int, minimum number of alleles required for a locus to be included
        returns:
            cleaned_data: pandas dataframe with the same columns as data with only the loci we wish to keep
            filtered_out: pandas dataframe with the loci that were filtered out so we can subset when loading the data
        '''
        # get the total allele counts for each locus
        locus_counts = cleaned_AC.groupby('LocusID').AC.sum()
        # filter out loci with less than min_alleles alleles
        filtered_out = locus_counts[locus_counts < min_alleles]

        # print('Number of loci filtered out:', len(filtered_out))

        cleaned_data = cleaned_AC[~cleaned_AC['LocusID'].isin(filtered_out.index)]

        # plot a histogram of the allele counts and plot the mean and median and print them in the legend rounded to the nearest integer
        plt.hist(cleaned_data.groupby("LocusID")['AC'].sum(), bins=100)
        plt.axvline(cleaned_data.groupby("LocusID")['AC'].sum().mean(), color='r', label='mean: '+str(round(cleaned_data.groupby("LocusID")['AC'].sum().mean())))
        plt.axvline(cleaned_data.groupby("LocusID")['AC'].sum().median(), color='g', label='median' + str(round(cleaned_data.groupby("LocusID")['AC'].sum().median())))
        plt.legend()
        plt.title('Histogram of total allele counts per locus, cleaned data')
        plt.xlabel('Total allele counts')
        plt.ylabel('Frequency')
        plt.show()

        return cleaned_data, filtered_out.index.tolist()
        
    @staticmethod
    def calculate_repeatmasker_span(data): # we still need this NOTE the locusid and allele_numbers are in terms of the AP tdb, need to convert first
        '''
        calculate the percentage of each allele that is made up of RepeatMasker annotations.

        Args: data (DataFrame): DataFrame with columns 'LocusID', 'allele_number', 'query_start', 'query_end', 'allele_length'
        Returns: DataFrame with an additional column 'allele_percent_element'
        '''
        # Discard rows with NaN in query_start and query_end
        data_copy = data.copy()
        data = data.dropna(subset=['query_start', 'query_end'])

        # Initialize a list to store the results
        results = []

        # Group by LocusID and allele_number to process each allele separately
        grouped = data.groupby(['LocusID', 'allele_number'])

        for name, group in grouped:
            # Sort the group by query_start
            group = group.sort_values(by='query_start')

            # Initialize variables to track the total covered length and the current annotation span
            total_covered_length = 0
            current_start = None
            current_end = None

            for _, row in group.iterrows():
                start = row['query_start']
                end = row['query_end']

                if current_start is None:
                    # Initialize the first annotation span
                    current_start = start
                    current_end = end
                else:
                    if start <= current_end:
                        # Overlapping annotation, extend the current span
                        current_end = max(current_end, end)
                    else:
                        # Non-overlapping annotation, add the current span to the total covered length
                        total_covered_length += current_end - current_start + 1
                        # Start a new span
                        current_start = start
                        current_end = end

            # Add the last span to the total covered length
            if current_start is not None:
                total_covered_length += current_end - current_start + 1

            # Calculate the percentage of the allele covered by RepeatMasker annotations
            allele_length = group['allele_length'].iloc[0]
            allele_percent_element = total_covered_length / allele_length

            # Append the result to the list
            results.append({
                'LocusID': name[0],
                'allele_number': name[1],
                'allele_percent_element': float(allele_percent_element)  # Ensure the value is a standard Python float
            })

        # Convert the results list to a DataFrame
        results_df = pd.DataFrame(results)

        # Merge the results back into the original DataFrame
        data = pd.merge(data_copy, results_df, on=['LocusID', 'allele_number'], how='outer')

        # Fill NaN values in allele_percent_element with 0
        data['allele_percent_element'].fillna(0, inplace=True)

        return data
    
    @staticmethod
    def harmonize_filters(data,allele_counts, new_afilters, new_lfilters): # TODO update calls to this, changed argument order
        '''
        Harmonize the allele and loci filters such that all alleles and loci that are not present in 
        AC are included in afilters and lfilters in addition to the new filters

        Args:
            allele_counts: updated allele counts that should reflect initial sfilters and allele filters
            new_afilters: dataframe of alleles that should be filtered according to allele_QC
            new_lfilters: list of loci that should be filtered according to filter_low_support_loci
            data: the data initially loaded with tdb.load_tdb that hasn't yet been filtered
        
        Returns:
            final_afilters: 2 column dataframe of LocusID allele_number pairs that should be filtered according to allele_QC and allele_counts
            final_lfilters: updated list of loci that should be filtered according to filter_low_support_loci and allele_counts
        '''
        # get all loci and alleles
        final_loci = allele_counts["LocusID"].unique()
        final_alleles = allele_counts[["LocusID","allele_number"]]
        final_alleles["LocusID:allele_number"] = final_alleles["LocusID"].astype(str) + ":" + final_alleles["allele_number"].astype(str)
        final_alleles = final_alleles.drop_duplicates()

        # get all alleles that are in data but not in AC TODO We need to format the afilters to allow for easily dropping locusID:allele pairs from the dataframe
        data['allele']["LocusID:allele_number"] = data['allele']["LocusID"].astype(str) + ":" + data['allele']["allele_number"].astype(str)
        allele_counts["LocusID:allele_number"] = allele_counts["LocusID"].astype(str) + ":" + allele_counts["allele_number"].astype(str)

        alleles_not_in_AC = data['allele'][~data['allele']["LocusID:allele_number"].isin(allele_counts["LocusID:allele_number"])]

        # get all alleles that are not currently in afilters and not in AC but are in data, add these to final_afilters
        final_afilters = pd.concat([new_afilters[["LocusID","allele_number"]], alleles_not_in_AC[["LocusID","allele_number"]]]).drop_duplicates() # only include unique locusID:allele pairs

        # get all loci that are not currently in lfilters and not in AC but are in data (ie we also want to filter these because there are no records of them in AC)
        loci_not_in_AC = data['locus'][~data['locus']["LocusID"].isin(allele_counts["LocusID"])][["LocusID"]].LocusID.tolist()
        # get all loci that are not currently in lfilters and not in AC but are in data, add these to final_lfilters
        final_lfilters = list(set(new_lfilters + loci_not_in_AC))

        return final_afilters, final_lfilters
        
    @staticmethod
    def run_all_QC(data_path, period_loci_path = "", period_stats_path="", repeatmasker_path="", initial_lfilters = [],min_spanning_reads = 10, min_allele_support = None,save_filters = False, filter_path = "./",existing_filters=False):
        '''
        TODO update this to reflect the new QC pipeline
        Load data, get all filters for alleles and loci, and apply them 

        parameters:
            data_path: path to the tdb data
            period_loci_path: path to the period loci file
            period_stats_path: path to the period stats file
            repeatmasker_path: path to the repeatmasker file
            initial_lfilters: list containing the initial lfilter field in this format [("LocusID", "in", [1,2,3])]
            min_spanning_reads: int, minimum number of spanning reads required for an allele to be included per sample
            min_allele_support: int, minimum number of alleles required for a locus to be included after QC
            save_filters: bool, whether to save the filters to a file
            filter_path: path to save the filters with output prefix
            existing_filters: tuple of existing filters to apply to the data if this has previously been ran
        returns:
            cleaned_data: pandas dataframe with the cleaned data
            final_afilters: 2 column dataframe of LocusID allele_number pairs that should be filtered according to allele_QC and allele_counts
            final_lfilters: updated list of loci that should be filtered according to filter_low_support_loci and allele_counts
        '''
        if not existing_filters:
            # Load data with minimum spanning reads and for the loci we are interested in (unfiltered)
            all_data = tdb.load_tdb(data_path, 
                            lfilters=initial_lfilters,sfilters=[('spanning_reads', '>=', min_spanning_reads)])
            print("all_data has loaded! Dimensions are: ", all_data["locus"].shape)
            
            # get allele counts and frequency for all samples so we can use them to filter loci
            allele_counts = TDBUtils.allele_count(all_data)

            # load additional data files for filtering steps
            period_stats = tdb_qc.load_period_data(period_loci_path, period_stats_path)
            
            # load repeatmasker data and generate the repeatmasker span column
            repeatmasker = tdb_qc.load_repeatmasker_data(repeatmasker_path)

            # combine the input data with the period, repeatmasker data, and allele data
            all_data_alleles = tdb_qc.combine_input_data(all_data['allele'], period_stats, repeatmasker)

            # calculate the percentage of each allele that is made up of RepeatMasker annotations
            all_data_alleles = tdb_qc.calculate_repeatmasker_span(all_data_alleles)

            # drop duplicates going forward
            all_data_alleles = all_data_alleles.drop_duplicates(subset=["LocusID","allele_number"])

            # perform allele QC based on the above stats
            # cleaned_data, new_afilters = allele_QC(all_data_alleles) # both of these returns are subsets of all_data_alleles
            res = tdb_qc.allele_QC(all_data_alleles)
            cleaned_data = res[0]
            new_afilters = res[1]

            # update allele counts after filtering
            cleaned_AC = tdb_qc.update_AC(allele_counts, new_afilters)

            # filter loci with low support
            cleaned_data, new_lfilters = tdb_qc.filter_low_support_loci(cleaned_AC, min_allele_support) # new_lfilters is a list of loci that were filtered out

            # harmonize the filters such that all alleles that are not present in AC are included in afilters and loci not in allele_counts are included in lfilters
            final_afilters, final_lfilters = tdb_qc.harmonize_filters(all_data,cleaned_AC, new_afilters, new_lfilters)

            if save_filters:
                print("Saving filters to ", filter_path + "afilters.tsv and ", filter_path + "lfilters.tsv")
                # save the afilters to a file at filter_path
                final_afilters.to_csv(filter_path + "afilters.tsv", index=False, sep="\t", header=False)
                # save the lfilters to a file at filter_path
                pd.Series(final_lfilters).to_csv(filter_path + "lfilters.tsv", index=False, header=False, sep="\t")
            # reload the data with the new lfilter -- this will filter loci across all tables
            del all_data # save RAM, get rid of this first

        else: # TODO add exception handeling for if the files don't exist
            # there are existing filters to load from files
            print("Loading existing filters from ", filter_path + "afilters.tsv and ", filter_path + "lfilters.tsv")
            final_afilters = pd.read_csv(filter_path + "afilters.tsv", sep="\t", header=None)
            final_afilters.columns = ["LocusID","allele_number"]
            final_lfilters = pd.read_csv(filter_path + "lfilters.tsv", sep="\t",header=None)[0].tolist()

        final_afilters["LocusID:allele_number"] = final_afilters["LocusID"].astype(str) + ":" + final_afilters["allele_number"].astype(str)
        final_lfilters_statement = [("LocusID", "not in", final_lfilters)]
        lfilters = initial_lfilters + final_lfilters_statement
        print("Loading data with final lfilters: ", lfilters)
        cleaned_data = tdb.load_tdb(data_path, lfilters=lfilters, sfilters=[('spanning_reads', '>=', min_spanning_reads)])

        # filter alleles that are in final_afilters from the allele table AND the sample table -- this is only for alleles that should be dropped but whose locus should be kept
        # this is necessary because there is not a good way to filter these when reading the tables in ie tdb doesn't allow for filtering based on 2 columns of the table

        # drop LocusID allele_number pairs from the allele table based on final_afilters dataframe
        print("Before filtering, cleaned_data has ", cleaned_data['allele'].shape[0], " alleles")
        cleaned_data['allele']['LocusID:allele_number'] = cleaned_data['allele']["LocusID"].astype(str) + ":" + cleaned_data['allele']["allele_number"].astype(str)
        cleaned_data['allele'] = cleaned_data['allele'][~cleaned_data['allele']["LocusID:allele_number"].isin(final_afilters["LocusID:allele_number"])]
        print("After filtering, cleaned_data has ", cleaned_data['allele'].shape[0], " alleles")
        # drop LocusID allele_number pairs from the sample table based on final_afilters dataframe
        for sample in cleaned_data['sample']:
            cleaned_data['sample'][sample]['LocusID:allele_number'] = cleaned_data['sample'][sample]["LocusID"].astype(str) + ":" + cleaned_data['sample'][sample]["allele_number"].astype(str)
            cleaned_data['sample'][sample] = cleaned_data['sample'][sample][~cleaned_data['sample'][sample]["LocusID:allele_number"].isin(final_afilters["LocusID:allele_number"])]

        return cleaned_data, final_afilters, final_lfilters

class tdb_MS_qc:
    '''
    This class is revised to be used when the MS data is included in the tdb object and can be used
    This is the class that MUST be used if using VC data because you need the motif spans to get
    indpendent allele lengths and sequences
    '''

    @staticmethod
    def log_changes(func):
        @functools.wraps(func)
        def wrapper_log_changes(*args, **kwargs):
            df = args[0] 
            
            # check if LocusID:allele_number is in columns
            if 'LocusID:allele_number' not in df.columns and "LocusID" in df.columns and "allele_number" in df.columns:
                df['LocusID:allele_number'] = df['LocusID'].astype(str) + ":" + df['allele_number'].astype(str)

            initial_alleles = df['LocusID:allele_number'].nunique() if 'LocusID:allele_number' in df.columns else 0
            initial_loci = df['LocusID'].nunique() if 'LocusID' in df.columns else 0
            
            result = func(*args, **kwargs)

            if isinstance(result, tuple):
                result_df = result[0]
            else:
                result_df = result
            
            # check if LocusID:allele_number is in columns
            if 'LocusID:allele_number' not in result_df.columns and "LocusID" in result_df.columns and "allele_number" in result_df.columns:
                result_df['LocusID:allele_number'] = result_df['LocusID'].astype(str) + ":" + result_df['allele_number'].astype(str)

            final_alleles = result_df['LocusID:allele_number'].nunique() if 'allele_number' in result_df.columns else 0
            final_loci = result_df['LocusID'].nunique() if 'LocusID' in result_df.columns else 0
            
            print(f"Function {func.__name__}:")
            print(f"  Alleles: {initial_alleles} -> {final_alleles}")
            print(f"  Loci: {initial_loci} -> {final_loci}")
            
            return result
        return wrapper_log_changes
    
    # Methods for processing and parsing MS data
    @staticmethod
    def add_purity_MS(data, *args, **kwargs): # this may suffice for both purity and MS, could probably also add methylation this way
        """
        Add motif span per allele to allele table from sample tables
        """
        allele = data['allele'].set_index(["LocusID", "allele_number"])
        parts = []
        for samp in data['sample']:
            parts.append((data['sample'][samp]
                        .set_index(["LocusID", "allele_number"])
                        .join(allele)
                        .reset_index()
                        [["LocusID", "allele_number", "allele_length",
                            "purity", "motif_span"]]
                        ))

        return pd.concat(parts).drop_duplicates()
    @staticmethod
    def get_tdb_mappings(tdb1, tdb2):
        '''
        Get the mappings between two tdbs to get locusid allele numbers
        This is so we can use any of the files we produced from a  previous tdb on the same regions
        but that have reordered allele numbers 
        '''
        check_locusIDs = pd.merge(tdb1["locus"], tdb2["locus"], on="LocusID", suffixes=("_AP", "_HPRC"), how="outer", indicator=True) 
        print("Checking the LocusIDs are the same for the two tdbs....")
        print(check_locusIDs[(check_locusIDs.start_AP != check_locusIDs.start_HPRC) | (check_locusIDs.end_AP != check_locusIDs.end_HPRC) | (check_locusIDs.chrom_AP != check_locusIDs.chrom_HPRC)]) # this says they should be the same, are the allele numbers different?
        print()
        # check the allele numbers
        check_alleles = pd.merge(tdb1["allele"], tdb2["allele"], on=["LocusID","sequence","allele_length"], suffixes=("_AP", "_HPRC"), how="outer", indicator=True)
        print("If all alleles in tdb1 are in tdb2, the merge values should all be both....")
        print(check_alleles._merge.value_counts()) # this should be 0 if they are the same
        #lets consolidate these in the repeatmasker data so we can map to the new numbers
        tdb_allele_mapping = pd.DataFrame({"LocusID":check_alleles.LocusID, "allele_number_MS":check_alleles.allele_number_HPRC, "allele_number_AP":check_alleles.allele_number_AP})
        return tdb_allele_mapping

    @staticmethod
    def split_sequence_MS(LocusID, allele_number, allele_len, sequence, motif_span, motif_list):
        """
        Split a sequence into the motif span
        """
        # for every motif span, split the sequence by the indeces designated in each (x-y) interval and relabel the sequence as LocuID:allele_number:motif:position
        split_sequences = {"LocusID":[],"allele_number":[], "sequence": [], "allele_len": [], "motif":[], "motif_start":[], "motif_end":[],"new_sequence":[],"new_allele_len":[]}
        # split the motif span into its intervals by "_"
        if motif_span == ".":
            # return the full sequence if no motif span is present
            split_sequences["LocusID"].append(LocusID)
            split_sequences["allele_number"].append(allele_number)
            split_sequences["motif"].append(motif_list)
            split_sequences["motif_start"].append(0)
            split_sequences["motif_end"].append(len(sequence))
            split_sequences["sequence"].append(sequence)
            split_sequences["allele_len"].append(allele_len)
            split_sequences["new_sequence"].append(sequence)
            split_sequences["new_allele_len"].append(allele_len)
            return pd.DataFrame(split_sequences)
        motif_spans = motif_span.split("_") if "_" in motif_span else [motif_span]
        # for each motif span in motif_spans, the format is motif_index(span_start, span_end), split the sequence by the span_start and span_end and key is LocusID:allele_number:motif_index:interval
        for i, span in enumerate(motif_spans):
            # get the start and end of the span and motif index
            motif_index, interval = span.split("(")
            motif_index = int(motif_index)
            start, end = interval.split("-")
            end = end[:-1]
            start = int(start)
            end = int(end) # inclusive
            # get the sequence for the span
            split_sequence = sequence[start:end+1]
            split_sequences["LocusID"].append(LocusID)
            split_sequences["allele_number"].append(allele_number)
            split_sequences["motif"].append(motif_list[motif_index])
            split_sequences["motif_start"].append(start)
            split_sequences["motif_end"].append(end)
            split_sequences["sequence"].append(sequence)
            split_sequences["allele_len"].append(allele_len)
            split_sequences["new_sequence"].append(split_sequence)
            split_sequences["new_allele_len"].append(len(split_sequence))

        return pd.DataFrame(split_sequences)
    
    @staticmethod
    def get_percent_spanned(allele_len, MS):
        # process the MS to get the total base pairs spanned
        # split the MS by "_"
        if MS == ".":
            return 0
        MS = str(MS).split("_") if "_" in str(MS) else [MS]
        # for each motif span in MS, the format is motif_index(span_start, span_end)
        # get the length of all ranges as span_end - span_start and sum them to get total base pairs spanned
        total_spanned = 0
        for span in MS:

            # get the start and end of the span and motif index
            motif_index, interval = span.split("(")
            start, end = interval.split("-")
            end = end[:-1]
            start = int(start)
            end = int(end)
            total_spanned += (end - start)
        # get the percent spanned
        percent_spanned = (total_spanned +1) / allele_len
        return percent_spanned

    @staticmethod
    def get_span_gaps(split):
        '''
        Get the length of the interruptions between motif spans
        Note: only input alleles that have multiple rows -- we should test with and without
        '''
        split["gap"] = 0
        # since we have sorted the data by LocusID, allele number, motif start and end we can just take the difference between the start of the next row and the end of the current row
        for LocusID, allele_number in split[["LocusID","allele_number"]].drop_duplicates().values:
            # get the rows for this LocusID and allele number
            rows = split[(split["LocusID"] == LocusID) & (split["allele_number"] == allele_number)]
            # get the pairwise differences between the start of the next row and the end of the current row
            pairwise_diffs = (rows["motif_start"].shift(-1) - rows["motif_end"]).iloc[:-1].tolist()
            split.loc[rows.index[:-1], "gap"] = pairwise_diffs
            split.loc[rows.index[-1], "gap"] = pairwise_diffs[-1] # set the last gap to the last value
        return split
    
    @staticmethod
    def _has_gaps(intervals):
        """
        Check if a list of intervals has gaps.
        Each interval is represented as a tuple (start, end).
        """
        # Sort intervals by start position
        intervals = sorted(intervals, key=lambda x: x[0])
        
        # Check for gaps
        for i in range(1, len(intervals)):
            if intervals[i][0] > intervals[i-1][1]:  # Gap exists
                return True
        return False

    @staticmethod
    def simple_span_group(MS):
        '''
        Straight forward span group assignment based only on MS
        Does not calculate gap, is allele specific instead of span specific
        should do the job for separating multispan from interrupted
        when we don't care about the gaps/ only care about including trimmed and multispan
        '''
        # (1) fully spanned by MS is designated previously using check_needs_split
        # (2) trimmed by MS, only one MS but new sequence is shorter than original sequence
        # (3) large interrupted MS, multiple MS with gaps and (4) multiple spans with no gaps
        
        # (1) check trimmed (only 2 parantheses in the string)
        if MS == ".":
            return "no_span"
        if MS.count("(") == 1:
            return "trimmed"
        # (2) to check interrupted vs multiple spans, we need to check if the spans have gaps or not
        # a gap is definted as a difference between the end of one span and the start of the next span
        # convert MS to a list of intervals -- the intervals are defined as (start-end) pairs separated by _i
        # where i is the index of the motif (we only care about the start and end of the interval)
        intervals = []
        for span in MS.split("_"):
            # get the start and end of the span and motif index
            print(span.split("("))
            motif_index, interval = span.split("(")
            start, end = interval.split("-")
            end = end[:-1]
            start = int(start)
            end = int(end) # inclusive
            intervals.append((start, end))
        # check if there are gaps
        if tdb_MS_qc._has_gaps(intervals):
            return "interrupted"
        else:
            return "multiple_spans"

    #TODO function to check if splitting the sequences corrected inclusion of element found by repeatmasker by checking the intervals in the sequence with motifs vs the intervals with the repeatmasker element
    @staticmethod
    def check_span_correction(split, repeatmasker):
        """
        Check if the span correction worked by checking the intervals in the sequence with motifs vs the intervals with the repeatmasker element
        """
        # get the repeatmasker intervals
        repeatmasker_intervals = repeatmasker[["LocusID", "allele_number", "query_start", "query_end"]].drop_duplicates()
        # get the split intervals
        split_intervals = split[["LocusID", "allele_number","allele_len", "motif_start", "motif_end"]].drop_duplicates()
        # merge the two dataframes on LocusID and allele number
        merged = pd.merge(repeatmasker_intervals, split_intervals, on=["LocusID", "allele_number"])
        # check if the intervals overlap
        merged["overlap"] = (merged["query_start"] <= merged["motif_end"]) & (merged["query_end"] >= merged["motif_start"])
        merged["original_overlap_percent"] = (merged["query_end"] - merged["query_start"])/merged.allele_len
        merged["new_overlap_percent"] = (merged[["query_end", "motif_end"]].min(axis=1) - merged[["query_start", "motif_start"]].max(axis=1)).clip(lower=0) / (merged["motif_end"] - merged["motif_start"])

        # plot the correction
        plt.hist(merged["original_overlap_percent"], bins=100, alpha=0.5, label="Original Overlap Percent")
        plt.hist(merged["new_overlap_percent"], bins=100, alpha=0.5, label="New Overlap Percent")
        plt.xlabel("Overlap Percent")
        plt.ylabel("Count")
        plt.legend()
        plt.title("Overlap Percent Before and After Span Correction")
        plt.show()
        return merged
    #TODO function that separates and labels alleles as being (1) fully spanned by MS (2) trimmed by MS (single MS) (3) large interrupted MS (4) multiple motifs without large interruptions
    @staticmethod
    def get_span_groups(split, interruption_threshold=0):
        # (1) fully spanned by MS, sequence is equal to new sequence
        fully_spanned = split[split["new_sequence"] == split["sequence"]].copy()
        fully_spanned["span_group"] = "fully_spanned"

        diff_seqs = split[split["new_sequence"] != split["sequence"]]
        multi_MS = diff_seqs[diff_seqs[["LocusID","allele_number"]].duplicated(keep=False)]
        trimmed = diff_seqs[~diff_seqs[["LocusID","allele_number"]].duplicated(keep=False)]
        # (2) trimmed by MS, only one MS but new sequence is shorter than original sequence
        trimmed["span_group"] = "trimmed"
        # (3) large interrupted MS, multiple MS with gaps and (4) multiple spans with no gaps
        # get gaps for these groups to separate into interrupted and multiple motifs
        multi_MS = multi_MS.sort_values(by=["LocusID","allele_number","motif_start","motif_end"])
        multi_MS = multi_MS.reset_index(drop=True)
        multi_MS = tdb_MS_qc.get_span_gaps(multi_MS)
        # get large interruptions -- TODO define a threshold for large interruptions
        interrupted = multi_MS[multi_MS["gap"] > interruption_threshold]
        interrupted["span_group"] = "interrupted"
        # get multiple spans with no gaps
        multi_MS = multi_MS[multi_MS["gap"] <= interruption_threshold]
        multi_MS["span_group"] = "multiple_spans"

        # combine all the groups
        all_groups = pd.concat([fully_spanned, trimmed, interrupted, multi_MS])
        # fill gap with 0
        all_groups["gap"] = all_groups["gap"].fillna(0)
        return all_groups
    
    @staticmethod
    def parse_vcid(vcid):
        # the VC is defined by an ID which is a list of TRIDS starting with ID= then delimited by ',' followed by ;MOTIFS= then delimited by ',' followed by ;STRUC
        # first split by ;
        vcid = vcid.split(';')
        # get the ID
        ID = vcid[0].split('=')[1]
        # split the ID by ',' to get individual TRIDs
        ID = ID.split(',')
        # get the motifs
        MOTIFS = vcid[1].split('=')[1]
        # split the motifs by ',' to get individual motifs
        MOTIFS = MOTIFS.split(',')
        # we don't really care about the structure
        return ID
    @staticmethod
    def parse_vcid_motifs(vcid):
        # the VC is defined by an ID which is a list of TRIDS starting with ID= then delimited by ',' followed by ;MOTIFS= then delimited by ',' followed by ;STRUC
        # first split by ;
        vcid = vcid.split(';')
        # get the ID
        ID = vcid[0].split('=')[1]
        # split the ID by ',' to get individual TRIDs
        ID = ID.split(',')
        # get the motifs
        MOTIFS = vcid[1].split('=')[1]
        # split the motifs by ',' to get individual motifs
        MOTIFS = MOTIFS.split(',')
        # we don't really care about the structure
        return MOTIFS
    @staticmethod
    def get_motif_periods(motifs):
        # get the period of each motif
        periods = []
        for motif in motifs:
            period = len(motif)
            periods.append(period)
        return periods

    @staticmethod
    def load_period_data(period_loci_path, period_stats_path): # TODO this will be replaced with MS, won't need external data
        '''
        Load the period data and return a pandas dataframe 
        '''
        period_locusids = pd.read_csv(period_loci_path, sep=' ', header=None)
        period_locusids.columns = ["period", "LocusID"]
        period_stats = pd.read_csv(period_stats_path, sep="\t")
        period_stats.columns = ["LocusID:allele", "longest_period","shortest_period","most_abundant_period","percent_input_match"]
        #split the LocusID:allele column into LocusID and allele columns
        period_stats[["LocusID","allele_number"]] = period_stats["LocusID:allele"].str.split(":", expand=True)
        period_stats["LocusID"] = period_stats['LocusID'].astype(int)
        period_stats["allele_number"] = period_stats['allele_number'].astype(int)

        # merge period_stats with period_locusids
        period_stats = pd.merge(period_stats, period_locusids, on="LocusID", how = "left")

        # print("There are ", period_stats.shape[0], " alleles in the period_stats dataframe")
        # print("There are ", period_stats["LocusID"].nunique(), " loci in the period_stats dataframe")

        return period_stats
    
    @staticmethod
    def load_repeatmasker_data(repeatmasker_path):
        '''
        Load the repeatmasker data and return a pandas dataframe
        '''
        # repeatmasker inputs -- NOTE there are NOT repeatmasker outputs for every locus in the database ie entries that don't have repeatmasker annotations will not be included
        repeatmasker = pd.read_csv(repeatmasker_path, sep="\t", header=None)
        # drop the first row
        repeatmasker = repeatmasker.dropna()
        repeatmasker.columns = ["LocusID:allele","query_start","query_end","subfamily","family"]
        # split the LocusID:allele column
        repeatmasker["LocusID"] = repeatmasker["LocusID:allele"].str.split(":", expand=True)[0].astype(int)
        repeatmasker["allele_number"] = repeatmasker["LocusID:allele"].str.split(":", expand=True)[1].astype(int)

        # calulate the length of the element detected by repeatmasker calculate_repeatmasker_span FIXME I dont think this is necessary
        repeatmasker["element_length"] = repeatmasker["query_end"].astype(int) - repeatmasker["query_start"].astype(int)


        return repeatmasker

    @staticmethod
    @log_changes
    def combine_input_data(data, period_data, repeatmasker_data): # TODO update to only include the repeatmasker data, may not need at all
        '''
        Combine the input data with the period and repeatmasker data
        parameters:
            data: data['allele'] table, pandas dataframe with columns 'LocusID','allele_number', 'percent_input_match','percent_element', 'allele_length', and 'period'
            period_data: pandas dataframe with columns 'LocusID','allele_number', 'longest_period','shortest_period','most_abundant_period','percent_input_match'
            repeatmasker_data: pandas dataframe with columns 'LocusID','allele_number','query_start','query_end','subfamily','family','element_length'
        '''
        # merge the period data with the input data
        data = pd.merge(data, period_data, on=["LocusID","allele_number"], how = "left")
        data["LocusID:allele"] = data["LocusID"].astype(str) + ":" + data["allele_number"].astype(str)
        # print("all_data after merging with period_stats data has  ", data["LocusID:allele"].nunique(), " alleles")
        # print("all_data after merging with period_stats data has  ", data["LocusID"].nunique(), " loci")
        # merge the repeatmasker data with the input data
        data = pd.merge(data, repeatmasker_data, on=["LocusID","allele_number"], how = "left")
        data["LocusID:allele"] = data["LocusID"].astype(str) + ":" + data["allele_number"].astype(str)
        # print("all_data after merging with repeatmasker data has  ", data["LocusID:allele"].nunique(), " alleles")
        # print("all_data after merging with repeatmasker data has  ", data["LocusID"].nunique(), " loci")

        return data

    @staticmethod
    @log_changes
    def allele_QC(data, exclude_allele_QC = []): # TODO update to include exclusions -- this should work fine with updated percent_match inputs but for now just allow exclusions
        '''
        Filter out alleles based on their sequence content including non-TR region inclusion
        parameters:
            data: pandas dataframe with columns 'LocusID','allele_number', 'percent_input_match','percent_element', 'allele_length', and 'period'
            exclude_allele_QC: list of LocusIDs to exclude from the QC
        returns:
            cleaned_data: pandas dataframe with the same columns as data 
            to_drop: pandas dataframe with the alleles that were filtered out so we can subset when loading the data
        '''
        # filter 1: remove alleles with high percent repeatmasker annotation and low percent input match (TE inclusion)
        data_2 = data[np.where((data["allele_percent_element"] <=0) | (data["percent_input_match"] > 60) | (data["LocusID"].isin(exclude_allele_QC)),True, False)]
        # filter 2: remove alleles that are outliers in terms of percent input period (neighboring TR and flanking sequence inclusion)

        #fill -inf values with 0
        data_2['percent_input_match'] = data_2['percent_input_match'].replace(-np.inf, 0)
        # calculate the outliers using IQR
        Q1 = data_2['percent_input_match'].quantile(0.25)
        Q3 = data_2['percent_input_match'].quantile(0.75)
        IQR = Q3 - Q1

        outliers_percent_input_match = data_2[(data_2['percent_input_match'] < (Q1 - 1.5 * IQR))] # only filter out the lower outliers -- should allow 100% input match alleles to be included

        # plot the distribution of percent input match and a red line at the cutoff
        plt.hist(data_2['percent_input_match'], bins=100)
        plt.axvline(x=(Q1 - 1.5 * IQR), color='r', label='IQR cutoff')
        # add labels and title
        plt.xlabel('Percent input match')
        plt.ylabel('Frequency')
        plt.title('Distribution of percent input match')
        plt.legend()
        plt.show()

        # FIXME the final outliers should be based on ALL filtered alleles
        final_outliers = outliers_percent_input_match[~((outliers_percent_input_match.allele_length < 40)&(outliers_percent_input_match.allele_length > outliers_percent_input_match.period))&(~outliers_percent_input_match.LocusID.isin(exclude_allele_QC))] # filter these from the outlier list
        clean_data = data_2[~data_2.index.isin(final_outliers.index)] # filter these from the non-outlier list
        
        to_drop = data[~data["LocusID:allele"].isin(clean_data["LocusID:allele"])] # get the alleles that were filtered out
        
        # ensure that there are no alleles with locusid in exclude_allele_QC in the to_drop list
        assert len(to_drop[to_drop["LocusID"].isin(exclude_allele_QC)]) == 0, "There are alleles in the to_drop list that are in the exclude_allele_QC list, something went wrong"
        
        return clean_data, to_drop

    @staticmethod
    @log_changes
    def update_AC(allele_counts, to_drop): # still good, will update what is in to_drop
        """
        Update allele counts after filtering by droppping alleles in to_drop

        Parameters:
        allele_counts: pandas DataFrame with columns 'LocusID', 'allele_number', 'AC', 'AF'
        to_drop: pandas DataFrame with columns 'LocusID', 'allele_number' returned by allele_QC

        Returns:
        clean_AC: cleaned allele_counts DataFrame
        """
        clean_AC = allele_counts[~allele_counts.set_index(['LocusID', 'allele_number']).index.isin(to_drop.set_index(['LocusID', 'allele_number']).index)]
        return clean_AC

    @staticmethod
    @log_changes
    def filter_low_support_loci(cleaned_AC, min_alleles):
        '''
        Filter out loci with less than min_alleles alleles in clean_AC and return the loci that are filtered out
        parameters:
            cleaned_AC: pandas dataframe with columns 'LocusID','allele_number', 'AC', and 'AF'
            min_alleles: int, minimum number of alleles required for a locus to be included
        returns:
            cleaned_data: pandas dataframe with the same columns as data with only the loci we wish to keep
            filtered_out: pandas dataframe with the loci that were filtered out so we can subset when loading the data
        '''
        # get the total allele counts for each locus
        locus_counts = cleaned_AC.groupby('LocusID').AC.sum()
        # filter out loci with less than min_alleles alleles
        filtered_out = locus_counts[locus_counts < min_alleles]

        # print('Number of loci filtered out:', len(filtered_out))

        cleaned_data = cleaned_AC[~cleaned_AC['LocusID'].isin(filtered_out.index)]

        # plot a histogram of the allele counts and plot the mean and median and print them in the legend rounded to the nearest integer
        plt.hist(cleaned_data.groupby("LocusID")['AC'].sum(), bins=100)
        plt.axvline(cleaned_data.groupby("LocusID")['AC'].sum().mean(), color='r', label='mean: '+str(round(cleaned_data.groupby("LocusID")['AC'].sum().mean())))
        plt.axvline(cleaned_data.groupby("LocusID")['AC'].sum().median(), color='g', label='median' + str(round(cleaned_data.groupby("LocusID")['AC'].sum().median())))
        plt.legend()
        plt.title('Histogram of total allele counts per locus, cleaned data')
        plt.xlabel('Total allele counts')
        plt.ylabel('Frequency')
        plt.show()

        return cleaned_data, filtered_out.index.tolist()
        
    @staticmethod
    def calculate_repeatmasker_span(data): # we still need this NOTE the locusid and allele_numbers are in terms of the AP tdb, need to convert first
        '''
        calculate the percentage of each allele that is made up of RepeatMasker annotations.

        Args: data (DataFrame): DataFrame with columns 'LocusID', 'allele_number', 'query_start', 'query_end', 'allele_length'
        Returns: DataFrame with an additional column 'allele_percent_element'
        '''
        # Discard rows with NaN in query_start and query_end
        data_copy = data.copy()
        data = data.dropna(subset=['query_start', 'query_end'])

        # Initialize a list to store the results
        results = []

        # Group by LocusID and allele_number to process each allele separately
        grouped = data.groupby(['LocusID', 'allele_number'])

        for name, group in grouped:
            # Sort the group by query_start
            group = group.sort_values(by='query_start')

            # Initialize variables to track the total covered length and the current annotation span
            total_covered_length = 0
            current_start = None
            current_end = None

            for _, row in group.iterrows():
                start = row['query_start']
                end = row['query_end']

                if current_start is None:
                    # Initialize the first annotation span
                    current_start = start
                    current_end = end
                else:
                    if start <= current_end:
                        # Overlapping annotation, extend the current span
                        current_end = max(current_end, end)
                    else:
                        # Non-overlapping annotation, add the current span to the total covered length
                        total_covered_length += current_end - current_start + 1
                        # Start a new span
                        current_start = start
                        current_end = end

            # Add the last span to the total covered length
            if current_start is not None:
                total_covered_length += current_end - current_start + 1

            # Calculate the percentage of the allele covered by RepeatMasker annotations
            allele_length = group['allele_length'].iloc[0]
            allele_percent_element = total_covered_length / allele_length

            # Append the result to the list
            results.append({
                'LocusID': name[0],
                'allele_number': name[1],
                'allele_percent_element': float(allele_percent_element)  # Ensure the value is a standard Python float
            })

        # Convert the results list to a DataFrame
        results_df = pd.DataFrame(results)

        # Merge the results back into the original DataFrame FIXME this may be where we are dropping our loci
        data = pd.merge(data_copy, results_df, on=['LocusID', 'allele_number'], how='outer')

        # Fill NaN values in allele_percent_element with 0
        data['allele_percent_element'].fillna(0, inplace=True)

        return data
    
    @staticmethod
    def harmonize_filters(data,allele_counts, new_afilters, new_lfilters): # TODO update calls to this, changed argument order
        '''
        Harmonize the allele and loci filters such that all alleles and loci that are not present in 
        AC are included in afilters and lfilters in addition to the new filters

        Args:
            allele_counts: updated allele counts that should reflect initial sfilters and allele filters
            new_afilters: dataframe of alleles that should be filtered according to allele_QC
            new_lfilters: list of loci that should be filtered according to filter_low_support_loci
            data: the data initially loaded with tdb.load_tdb that hasn't yet been filtered
        
        Returns:
            final_afilters: 2 column dataframe of LocusID allele_number pairs that should be filtered according to allele_QC and allele_counts
            final_lfilters: updated list of loci that should be filtered according to filter_low_support_loci and allele_counts
        '''
        # get all loci and alleles
        final_loci = allele_counts["LocusID"].unique()
        final_alleles = allele_counts[["LocusID","allele_number"]]
        final_alleles["LocusID:allele_number"] = final_alleles["LocusID"].astype(str) + ":" + final_alleles["allele_number"].astype(str)
        final_alleles = final_alleles.drop_duplicates()

        # get all alleles that are in data but not in AC TODO We need to format the afilters to allow for easily dropping locusID:allele pairs from the dataframe
        data['allele']["LocusID:allele_number"] = data['allele']["LocusID"].astype(str) + ":" + data['allele']["allele_number"].astype(str)
        allele_counts["LocusID:allele_number"] = allele_counts["LocusID"].astype(str) + ":" + allele_counts["allele_number"].astype(str)

        alleles_not_in_AC = data['allele'][~data['allele']["LocusID:allele_number"].isin(allele_counts["LocusID:allele_number"])]

        # get all alleles that are not currently in afilters and not in AC but are in data, add these to final_afilters
        final_afilters = pd.concat([new_afilters[["LocusID","allele_number"]], alleles_not_in_AC[["LocusID","allele_number"]]]).drop_duplicates() # only include unique locusID:allele pairs

        # get all loci that are not currently in lfilters and not in AC but are in data (ie we also want to filter these because there are no records of them in AC)
        loci_not_in_AC = data['locus'][~data['locus']["LocusID"].isin(allele_counts["LocusID"])][["LocusID"]].LocusID.tolist()
        # get all loci that are not currently in lfilters and not in AC but are in data, add these to final_lfilters
        final_lfilters = list(set(new_lfilters + loci_not_in_AC))

        return final_afilters, final_lfilters
    
    @staticmethod
    def combine_locus(new_tdb, curr_tdb):
        # drop locusids in curr_tdb that are in the new_tdb original locusids
        curr_tdb["locus"] = curr_tdb["locus"][~curr_tdb["locus"]["LocusID"].isin(new_tdb["locus"]["original_LocusID"])]
        # add original_locusid column to current tdb
        if "original_LocusID" not in curr_tdb["locus"].columns:
            curr_tdb["locus"]["original_LocusID"] = curr_tdb["locus"]["LocusID"]
        # concat the locus table from new_tdb to curr_tdb
        curr_tdb["locus"] = pd.concat([curr_tdb["locus"], new_tdb["locus"]], ignore_index=True).reset_index(drop=True)
        return curr_tdb
    @staticmethod
    def combine_allele(new_tdb, curr_tdb):
        # drop the allele rows with locusids in the new_tdb original locusids
        curr_tdb["allele"] = curr_tdb["allele"][~curr_tdb["allele"]["LocusID"].isin(new_tdb["locus"]["original_LocusID"])]
        # add original_locusid column to current tdb
        if "original_LocusID" not in curr_tdb["allele"].columns:
            curr_tdb["allele"]["original_LocusID"] = curr_tdb["allele"]["LocusID"]
        # concat the allele table from new_tdb to curr_tdb
        curr_tdb["allele"] = pd.concat([curr_tdb["allele"], new_tdb["allele"]], ignore_index=True).reset_index(drop=True)
        # expliticly cast all allele_number columns to int -- floating point are messing up the ids
        curr_tdb["allele"]["allele_number"] = curr_tdb["allele"]["allele_number"].astype(int)
        return curr_tdb
    @staticmethod
    def combine_sample(new_tdb, curr_tdb):
        # for each sample in the curr_tdb, drop the rows with locusids in the new_tdb original locusids and concat the new_tdb sample table
        for sample in curr_tdb["sample"].keys():
            curr_tdb["sample"][sample] = curr_tdb["sample"][sample][~curr_tdb["sample"][sample]["LocusID"].isin(new_tdb["locus"]["original_LocusID"])]
            # concat the sample table from new_tdb to curr_tdb
            curr_tdb["sample"][sample] = pd.concat([curr_tdb["sample"][sample], new_tdb["sample"][sample]], ignore_index=True).reset_index(drop=True)

            # set all allele_number columns to int -- floating point are messing up the ids
            curr_tdb["sample"][sample]["allele_number"] = curr_tdb["sample"][sample]["allele_number"].astype(int)


        return curr_tdb
    @staticmethod
    def append_tdb(tdb_pkl, curr_tdb):
        '''
        reads a pickle compressed tdb file produced by update_tables 
        and consolidates the data with the current tdb
        '''
        # load the tdb from pkl
        with open(tdb_pkl, "rb") as f:
            new_tdb = pd.read_pickle(f)
        curr_tdb = tdb_MS_qc.combine_locus(new_tdb, curr_tdb)
        # update the allele table
        curr_tdb = tdb_MS_qc.combine_allele(new_tdb, curr_tdb)
        # update the sample tables
        curr_tdb = tdb_MS_qc.combine_sample(new_tdb, curr_tdb)
        return curr_tdb
    
    @staticmethod
    def run_all_QC(data_path, period_loci_path = "", period_stats_path="", repeatmasker_path="", initial_lfilters = [],min_spanning_reads = 10,
                    min_allele_support = None,save_filters = False, filter_path = "./",existing_filters=False, update_tdb= None, exclude_allele_QC = []):
        '''
        Load data, get all filters for alleles and loci, and apply them 

        parameters:
            data_path: path to the tdb data
            period_loci_path: path to the period loci file
            period_stats_path: path to the period stats file
            repeatmasker_path: path to the repeatmasker file
            initial_lfilters: list containing the initial lfilter field in this format [("LocusID", "in", [1,2,3])]
            NOTE if using update_tdb, include the original ids in the initial_lfilters to INCLUDE if they aren't already
            min_spanning_reads: int, minimum number of spanning reads required for an allele to be included per sample
            min_allele_support: int, minimum number of alleles required for a locus to be included after QC
            save_filters: bool, whether to save the filters to a file
            filter_path: path to save the filters with output prefix
            existing_filters: tuple of existing filters to apply to the data if this has previously been ran
            update_tdb: path to pkl saved tdb dictionary to update the database with, if None nothing will be done
            exclude_allele_QC: list of LocusIDs to exclude from allele QC, if None nothing will be done differently

            NOTE exclude_allele_QC is meant for adding TRs that are from update_tdb that may not have entries in the percent_match data
            or whose entries are currently incorrect due to the allele sequences needing to be split/trimmed
            I might get rid of this if we have the correct data for the regions we want to include
            NOTE also the spannin reads filter may not actually be applied to the updated alleles so we may want to do that separately if we do 
            want to include it -- AC filters should still work on the updated alleles/ loci since it just counts how many times alleles come up 
            given the sample tables
        returns:
            cleaned_data: pandas dataframe with the cleaned data
            final_afilters: 2 column dataframe of LocusID allele_number pairs that should be filtered according to allele_QC and allele_counts
            final_lfilters: updated list of loci that should be filtered according to filter_low_support_loci and allele_counts
        '''
        if not existing_filters:
            # Load data with minimum spanning reads and for the loci we are interested in (unfiltered)
            all_data = tdb.load_tdb(data_path, 
                    lfilters=initial_lfilters,sfilters=[('spanning_reads', '>=', min_spanning_reads)])
            logging.info("Loaded all_data with dimensions: %s", all_data["locus"].shape)
            
            if update_tdb is not None:
                all_data = tdb_MS_qc.append_tdb(update_tdb, all_data) # this will update the database with the new data
                logging.info("Updated tdb with new data from %s", update_tdb)
                logging.info("Updated all_data dimensions: %s", all_data["locus"].shape)
                # head the all_data allele table to check the allele numbers are int
                print(all_data["allele"].head())
                # check how many of the excluded loci are in the new data
                if exclude_allele_QC is not None:
                    excluded_loci = all_data["locus"][all_data["locus"]["LocusID"].isin(exclude_allele_QC)]
                    logging.info("Excluded loci in new data: %s", excluded_loci.shape[0])

            # Get allele counts and frequency for all samples so we can use them to filter loci
            allele_counts = TDBUtils.allele_count(all_data)
            logging.info("Calculated allele counts")

            # Load additional data files for filtering steps
            period_stats = tdb_MS_qc.load_period_data(period_loci_path, period_stats_path)
            logging.info("Loaded period data")

            repeatmasker = tdb_MS_qc.load_repeatmasker_data(repeatmasker_path)
            logging.info("Loaded repeatmasker data")

            # Combine the input data with the period, repeatmasker data, and allele data
            all_data_alleles = tdb_MS_qc.combine_input_data(all_data['allele'], period_stats, repeatmasker)
            logging.info("Combined input data")
            if exclude_allele_QC is not None:
                    excluded_loci = all_data_alleles[all_data_alleles["LocusID"].isin(exclude_allele_QC)]
                    logging.info("Excluded loci in new data: %s", len(excluded_loci.LocusID.unique()))

            # Calculate the percentage of each allele that is made up of RepeatMasker annotations
            all_data_alleles = tdb_MS_qc.calculate_repeatmasker_span(all_data_alleles)
            logging.info("Calculated RepeatMasker span")
            if exclude_allele_QC is not None:
                    excluded_loci = all_data_alleles[all_data_alleles["LocusID"].isin(exclude_allele_QC)]
                    logging.info("Excluded loci in new data: %s", len(excluded_loci.LocusID.unique()))

            # Drop duplicates going forward
            all_data_alleles = all_data_alleles.drop_duplicates(subset=["LocusID","allele_number"])
            logging.info("Dropped duplicate alleles")
            if exclude_allele_QC is not None:
                    excluded_loci = all_data_alleles[all_data_alleles["LocusID"].isin(exclude_allele_QC)]
                    logging.info("Excluded loci in new data: %s", len(excluded_loci.LocusID.unique()))

            # Perform allele QC based on the above stats
            res = tdb_MS_qc.allele_QC(all_data_alleles, exclude_allele_QC)
            cleaned_data = res[0]
            new_afilters = res[1]
            logging.info("Performed allele QC")
            if exclude_allele_QC is not None:
                    excluded_loci = cleaned_data[cleaned_data["LocusID"].isin(exclude_allele_QC)]
                    logging.info("Excluded loci in new data: %s", len(excluded_loci.LocusID.unique()))

            # Update allele counts after filtering
            cleaned_AC = tdb_MS_qc.update_AC(allele_counts, new_afilters)
            logging.info("Updated allele counts after filtering")
            if exclude_allele_QC is not None:
                    excluded_loci = cleaned_AC[cleaned_AC["LocusID"].isin(exclude_allele_QC)]
                    logging.info("Excluded loci in new data: %s", excluded_loci.shape[0])

            # Filter loci with low support
            cleaned_data, new_lfilters = tdb_MS_qc.filter_low_support_loci(cleaned_AC, min_allele_support)
            logging.info("Filtered loci with low support")
            if exclude_allele_QC is not None:
                    excluded_loci = cleaned_data[cleaned_data["LocusID"].isin(exclude_allele_QC)]
                    logging.info("Excluded loci in new data: %s", len(excluded_loci.LocusID.unique()))

            # Harmonize the filters
            final_afilters, final_lfilters = tdb_MS_qc.harmonize_filters(all_data, cleaned_AC, new_afilters, new_lfilters)
            logging.info("Harmonized filters")

            if save_filters:
                logging.info("Saving filters to %safilters.tsv and %slfilters.tsv", filter_path, filter_path)
                final_afilters.to_csv(filter_path + "afilters.tsv", index=False, sep="\t", header=False)
                pd.Series(final_lfilters).to_csv(filter_path + "lfilters.tsv", index=False, header=False, sep="\t")

            # Reload the data with the new lfilter
            del all_data  # Save RAM
            logging.info("Deleted all_data to save RAM")

        else:
            logging.info("Loading existing filters from %safilters.tsv and %slfilters.tsv", filter_path, filter_path)
            final_afilters = pd.read_csv(filter_path + "afilters.tsv", sep="\t", header=None)
            final_afilters.columns = ["LocusID", "allele_number"]
            final_lfilters = pd.read_csv(filter_path + "lfilters.tsv", sep="\t", header=None)[0].tolist()

        final_afilters["LocusID:allele_number"] = final_afilters["LocusID"].astype(str) + ":" + final_afilters["allele_number"].astype(str)
        final_lfilters_statement = [("LocusID", "not in", final_lfilters)]
        lfilters = initial_lfilters + final_lfilters_statement
        logging.info("Loading data with final lfilters: %s", lfilters)
        cleaned_data = tdb.load_tdb(data_path, lfilters=lfilters, sfilters=[('spanning_reads', '>=', min_spanning_reads)])
        logging.info("Loaded cleaned data")
        
        # TODO append the update table again
        if update_tdb is not None:
                cleaned_data = tdb_MS_qc.append_tdb(update_tdb, cleaned_data) # this will update the database with the new data
                logging.info("Updated tdb with new data from %s", update_tdb)
                logging.info("Updated all_data dimensions: %s", cleaned_data["locus"].shape)
        # TODO drop loci from update that haven't been dropped
        cleaned_data["locus"] = cleaned_data["locus"][~cleaned_data["locus"]["LocusID"].isin(final_lfilters)]
        # TODO continue as normal, the alleles should be in afilters if they need to be filtered out
        # Filter alleles that are in final_afilters from the allele table and the sample table
        logging.info("Before filtering, cleaned_data has %d alleles", cleaned_data['allele'].shape[0])
        cleaned_data['allele']['LocusID:allele_number'] = cleaned_data['allele']["LocusID"].astype(str) + ":" + cleaned_data['allele']["allele_number"].astype(str)
        cleaned_data['allele'] = cleaned_data['allele'][~cleaned_data['allele']["LocusID:allele_number"].isin(final_afilters["LocusID:allele_number"])]
        logging.info("After filtering, cleaned_data has %d alleles", cleaned_data['allele'].shape[0])

        for sample in cleaned_data['sample']:
            cleaned_data['sample'][sample]['LocusID:allele_number'] = cleaned_data['sample'][sample]["LocusID"].astype(str) + ":" + cleaned_data['sample'][sample]["allele_number"].astype(str)
            cleaned_data['sample'][sample] = cleaned_data['sample'][sample][~cleaned_data['sample'][sample]["LocusID:allele_number"].isin(final_afilters["LocusID:allele_number"])]
        logging.info("Filtered alleles from sample tables")


        return cleaned_data, final_afilters, final_lfilters
    
    @staticmethod
    def update_tables(split_seqs, orig_data,max_locusid = 4542944):
        '''
        TODO figure out when we want to deal with interruptions -- right now assume that everything we want to keep is in the split_seqs and data
        Given a dataframe of alleles that have been split into multiple TR sequences,
        (1) reassign each part to an existing or new LocusID
        (2) add a row for each new locus to the locus table
        (3) add a column that has the linked Locusid and the motif
        (4) propogate changes to allele table: merge the allele table with the locus table to add the new locusIDs and original locus IDs
        (5) for each original locus, add a row for every new locus
        (6) populate the sequence column with the corresponding subsequence from the allele given in allele_number
        (7) propogate to sample tables: merge on original allele number and original LocusID
        (8) clean up allele and sample tables: collapse duplicate alleles to the smallest allele number
        '''
        # make a copy of the data to avoid modifying the original data
        data = tdb_MS_qc.copy_tdb(orig_data)
        data["locus"]["original_LocusID"] = data["locus"]["LocusID"]
        data["locus"]["span_number"] = 1 # the original entry will always represent the first span
        split_seqs["period"] = split_seqs["motif"].str.len() # add the period to the split_seqs dataframe
        split_seqs = split_seqs[(split_seqs["period"] >1) & (split_seqs["period"] < 7)] 

        # add span_number so we can align the new loci to which section of the allele sequence it came from
        all_rows = []
        # for every LocusID allele_number pair in split_seqs, add a value that is seq(number_of_rows) to the span_number column 
        for LocusID, allele_number in split_seqs[["LocusID","allele_number"]].drop_duplicates().values:
            # get the rows for this LocusID and allele number
            rows = split_seqs[(split_seqs["LocusID"] == LocusID) & (split_seqs["allele_number"] == allele_number)]
            # add a column for the span number
            rows["span_number"] = range(1, len(rows)+1)

            # FIXME -- this is hardcoded to account for SCA27B insertion in a single allele -- if we want this to be more permenant or apply this to new regions or samples we will need to generalize this
            # this is just not a naiive problem to solve and I need more time than I want to spend on this to solve it
            if LocusID == 3240174 and allele_number == 27:
                # switch the span numbers for these 2 rows
                rows["span_number"] = [2,1]
            
            all_rows.append(rows)
        split_seqs_2 = pd.concat(all_rows)

        allele_value_counts = split_seqs_2[["LocusID", "allele_number"]].value_counts().reset_index()
        allele_value_counts.columns = ["LocusID", "allele_number", "count"]

        # Merge the value counts back into split_seqs
        to_iterate = pd.merge(split_seqs_2, allele_value_counts, on=["LocusID", "allele_number"], how="left")
        to_iterate = to_iterate[to_iterate["count"] > 1][["LocusID", "allele_number","count"]].groupby(["LocusID"])['count'].max().reset_index() # this is a list of LocusIDs and allele numbers that have multiple rows
        curr_new_locus = max_locusid
        for LocusID, count in to_iterate.values:
            # get the rows for this LocusID and allele number
            rows = split_seqs_2[(split_seqs_2["LocusID"] == LocusID)]
            # get the new LocusID
            new_LocusID = curr_new_locus + 1
            curr_locus = data["locus"].loc[data["locus"]["LocusID"] == LocusID]
            # add a new row to the locus table for each new LocusID with the same coordinates and the current locusID as the original locusID
            # print(LocusID, count)
            for i in range(count-1):
                # print(LocusID)
                data["locus"] = pd.concat([data["locus"], pd.DataFrame({"LocusID": [new_LocusID], "original_LocusID": [LocusID], "chrom": [curr_locus.chrom.iloc[0]], "start": [curr_locus.start.iloc[0]], "end": [curr_locus.end.iloc[0]],"span_number": i+2})], ignore_index=True)
                new_LocusID += 1
            curr_new_locus = new_LocusID - 1
        # add period to the locus table from split_seqs_2 -- ADDED FOR PREPROCESSING -- this WILL throw an error about LocusID later
        data["locus"] = pd.merge(data["locus"], split_seqs_2[["LocusID","period","motif","span_number"]], left_on = ["original_LocusID","span_number"],right_on=["LocusID","span_number"], how="left")
        # rename LocusID from locus table back to LocusID and drop the LocusID_y column
        data["locus"].drop(columns=["LocusID_y"], inplace=True)
        data["locus"].rename(columns={"LocusID_x":"LocusID"}, inplace=True)

        data["locus"] = data["locus"].drop_duplicates(subset=["LocusID","original_LocusID","span_number"]).reset_index(drop=True)

        data["allele"]["original_LocusID"] = data["allele"]["LocusID"]
        data["allele"].drop(columns=["LocusID"], inplace=True)
        test = pd.merge(data["allele"], data["locus"][["LocusID", "original_LocusID","span_number"]], on="original_LocusID", how="left")
        
        print(test)
        test2 = pd.merge(test, split_seqs_2[["LocusID","allele_number","new_sequence","new_allele_len","span_number","period"]], left_on=["original_LocusID","allele_number","span_number"], right_on=["LocusID","allele_number","span_number"], how="left") # add new sequence
        test2.drop(columns=["LocusID_y"], inplace=True)
        test2.rename(columns={"LocusID_x":"LocusID"}, inplace=True)

        # update the allele sequence and length
        test2.drop(columns=['sequence','allele_length'], inplace=True)
        test2.rename(columns={"new_sequence":"sequence","new_allele_len":"allele_length"}, inplace=True)
        data["allele"] = test2.copy()
        del test2
        # update the allele numbers
        test3 = data["allele"].drop_duplicates(subset=["LocusID","sequence","allele_length"]) # remove duplicates after the sequences have been subset

        # renumber the alleles sequentially after 0 (0 is reference) by length
        for locus in test3["LocusID"].drop_duplicates():
            # get the rows for this locus
            rows = test3[test3["LocusID"] == locus]
            ref_allele = rows[rows["allele_number"] == 0].copy()
            ref_allele["new_allele_number"] = 0
            # remove the reference allele from the rows
            rows = rows[rows["allele_number"] != 0]
            # for every number after 0, sort the rows by length and assign the new allele number sequentially
            rows = rows.sort_values(by=["allele_length"], ascending=False)
            rows["new_allele_number"] = range(1, len(rows)+1)
            # add the reference allele back to the rows
            rows = pd.concat([ref_allele, rows])
            # add the rows to the test3 dataframe
            test3 = pd.concat([test3[test3["LocusID"] != locus], rows])

        test4 = pd.merge(data["allele"], test3[["LocusID","sequence","new_allele_number"]], on=["LocusID","sequence"], how="left")

        # update the sample tables with the new LocusIDs and alleles by merging on original LocusID and allele number
        for sample in data["sample"].keys():
            locus_counts = pd.merge(data["sample"][sample], test4[["LocusID","original_LocusID","allele_number","new_allele_number"]], left_on=["LocusID","allele_number"],
                right_on = ["original_LocusID","allele_number"], how="left").drop_duplicates(subset=["LocusID_y","new_allele_number"]).sort_values(by=["LocusID_y","new_allele_number"])[["LocusID_y"]].value_counts() 
            # drop the duplicates
            locus_counts = locus_counts.reset_index()
            test5 = pd.merge(data["sample"][sample], test4[["LocusID","original_LocusID","allele_number","new_allele_number"]], left_on=["LocusID","allele_number"],
                    right_on = ["original_LocusID","allele_number"], how="left").drop_duplicates(subset=["LocusID_y","new_allele_number"])
            # for every locus that has a single allele number, duplicate the row in test5 for that allele number and locus
            for locus, count in locus_counts[locus_counts["count"] == 1].values:
                # get the rows for this locus
                rows = test5[test5["LocusID_y"] == locus]
                # duplicate the rows
                rows = pd.concat([rows, rows])
                # add the rows to the test5 dataframe
                test5 = pd.concat([test5[test5["LocusID_y"] != locus], rows])

            # set the new_allele_number to the new allele number
            test5["allele_number"] = test5["new_allele_number"]
            # drop the new_allele_number column
            test5.drop(columns=["new_allele_number"], inplace=True)
            # reassign the LocusID to the LocusID_y column and drop the LocusID_y, LocusID_x, original_LocusID columns
            test5.rename(columns={"LocusID_y":"LocusID"}, inplace=True)
            test5.drop(columns=["LocusID_x","original_LocusID"], inplace=True)
            test5.reset_index(drop=True, inplace=True)
            data["sample"][sample] = test5.copy()
            del test5

        # clean up allele table since we don't need the original LocusID or the allele number
        test4.dropna(subset="sequence", inplace=True) # drop rows with no sequence -- these are not real alleles, never observed or artifacts of placeholders for spans not observed in other alleles
        test4.drop(columns=["allele_number"], inplace=True)
        # update the allele number
        test4.rename(columns={"new_allele_number":"allele_number"}, inplace=True)
        data["allele"] = test4.copy()
        data["allele"].drop_duplicates(subset=["LocusID","allele_number","sequence","allele_length"], inplace=True)
        del test4

        return data
    
    @staticmethod
    def copy_tdb(data):
        curr_copy = {}
        curr_copy["locus"] = data["locus"].copy()
        curr_copy["allele"] = data["allele"].copy()
        curr_copy["sample"] = {}
        for sample in data["sample"].keys():
            curr_copy["sample"][sample] = data["sample"][sample].copy()
        return curr_copy