import tdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sb
from scipy.stats import pearsonr
from scipy.stats import mannwhitneyu
from scipy.stats import kruskal
import matplotlib.colors as mcolors
import functools
from collections import defaultdict
import os
import pyranges as pr

class TDBUtils:
    """
    Utility class for handling TDB data and performing various analyses.

    This class provides static methods for:
    - Counting alleles and calculating allele frequencies.
    - Calculating sequence composition polymorphism scores.
    - Calculating length polymorphism scores.
    - Retrieving observed allele lengths.
    - Calculating purity for each allele.

    all methods are static and take a TDB data object as input.
    """

    @staticmethod
    def allele_count(data, samples=None, *args, **kwargs): # need this to filter loci with low total allele counts
        """
        Allele counts and frequency
        """
        samp_data = data["sample"]
        if samples is None:
            samples = samp_data.keys()
        all_alleles = pd.concat([samp_data[_][["LocusID", "allele_number"]]
                                for _ in samples])
        lcnts = all_alleles["LocusID"].value_counts()
        acnts = (all_alleles.groupby(["LocusID"])["allele_number"]
                .value_counts()
                .rename("AC")
                .reset_index(level=1))
        acnts['AF'] = acnts['AC'] / lcnts
        return (data['locus'].set_index("LocusID")
                .join(acnts)
                .fillna(0)
                .astype({'allele_number': np.uint16, 'AC': np.uint16})).reset_index()
    
    @staticmethod
    def composition_polymorphism_score(data, min_af=0.01, kmer_len=5, min_freq=5, *args, **kwargs):
        """
        Calculate loci's sequence composition as mean jaccard index
        """
        a_cnts = TDBUtils.allele_count(data).reset_index().set_index(
            ["LocusID", "allele_number"])
        a_cnts['sequence'] = data['allele'].set_index(
            ["LocusID", "allele_number"])['sequence']
        result = (a_cnts.reset_index()
                .where(lambda x: x["AF"] >= min_af)
                .groupby(['LocusID'])[["sequence", "AC"]]
                .apply(lambda x:
                        tdb.alleles_jaccard_dist(x["sequence"].values, x["AC"].values,
                                                kmer_len, min_freq)))
        result.name = "comp_poly_score"
        return pd.concat([data['locus'].set_index("LocusID"), result], axis=1)

    @staticmethod
    def length_polymorphism_score(data, min_af=0.01, *args, **kwargs):
        """
        Number of distinct alleles by length per 100 samples for each locus
        """
        a_cnts = TDBUtils.allele_count(data).reset_index().set_index(
            ["LocusID", "allele_number"])
        a_cnts['allele_length'] = data['allele'].set_index(
            ["LocusID", "allele_number"])['allele_length']
        len_cnts = (a_cnts.reset_index()
                    .where(lambda x: x["AF"] >= min_af)
                    .groupby(['LocusID'])['allele_length']
                    .nunique() / (len(data['sample']) / 100))
        len_cnts.name = 'len_poly_score'
        return pd.concat([data['locus'].set_index('LocusID'), len_cnts], axis=1)

    @staticmethod
    def get_observed_allele_lengths(data, samples):
        """
        Given a tdb and a list of sample names
        return a DataFrame of LocusID, allele_number, allele_length, sample_name
        """
        observed_allele_lengths = []
        for name in samples:
            sample_table = data["sample"][name]
            
            # Create an index of the alleles present in the sample
            idx_table = pd.MultiIndex.from_frame(sample_table[["LocusID", "allele_number"]])
            # Subset the alleles and pull their lengths
            alleles = data["allele"].set_index(["LocusID", "allele_number"])
            result = alleles.loc[idx_table, ["allele_length"]]
            result['sample_name'] = name
            observed_allele_lengths.append(result)
        # Sorting and indexing here speeds up calculating the distribution per-locus later
        return (pd.concat(observed_allele_lengths)
                    .reset_index()
                    .sort_values("LocusID")
                    .set_index("LocusID"))
    
    @staticmethod
    def purity(data, *args, **kwargs):
        """
        modify the methyl function to join the allele table with the sample data to get purity for each allele
        given that the purity is in the sample data
        """
        allele = data['allele'].set_index(["LocusID", "allele_number"])
        parts = []
        for samp in data['sample']:
            parts.append((data['sample'][samp]
                        .set_index(["LocusID", "allele_number"])
                        .join(allele)
                        .reset_index()
                        [["LocusID", "allele_number", "allele_length",
                            "purity"]]
                        ))

        return pd.concat(parts).drop_duplicates()
     
    # TODO move add_MS to this class instead of QC

# TODO define the following class such that we can import functions from the trgt paper (ie their purity script to recalculate purity per sequence)
# and add scripts from the TR Catalog preprint (provided as well as just made from what they say in their methods like kmer counting method)
class SequenceAttributes: #TODO add sequence analsysis functions here ie the kmer counting and purity calculation scripts from other sources
    """
    Utility class for handling sequence data and performing various analyses.

    This class provides static methods for:
    - Counting kmers in a sequence.
    - Calculating the reverse complement of a DNA sequence.
    - Rotating a string.
    - Determining the possible repeat units of a motif.
    - Getting the lexigraphically smallest repeat unit for a motif.

    all methods are static and take a sequence as input.
    """

    @staticmethod
    def rev_comp(seq):
        """
        Determine the reverse complement of a DNA sequence
        Arguments: 
            seq: sequence to reverse complement
        Returns:
            reverse complement of seq
        """ 
        comp_dict = {'A': 'T', 'T': 'A', 'C':'G', 'G':'C','N':'N'}
        new_seq = ''
        for base in seq[::-1]:
            new_seq = new_seq + comp_dict[base]
        return(new_seq)
    
    @staticmethod
    def rotate_string(seq):
        """
        Determine the list of all rotations of a string
        Arguments: 
            seq: str to rotate
        Returns:
            List of rotations of seq
        """ 
        ru_list = [seq]
        for i in range(len(seq) - 1):
            ru_list.append(ru_list[-1][1:]+ru_list[-1][0])
        return(list(set(ru_list)))
    
    @staticmethod
    def possible_ru(ru):
        """
        Determine the list of identicial repeat units through rotation and 
        reverse complementation

        Arguments:
            ru: repeat unit string
        Returns:
            List of possible repeat units
        """
        ru = ru.upper()
        ru_list = SequenceAttributes.rotate_string(ru)
        ru_list_rev_comp = []
        for item in ru_list:
            ru_list_rev_comp.append(SequenceAttributes.rev_comp(item))
        final_list = list(set(ru_list + ru_list_rev_comp))
        return(final_list)
    
    @staticmethod
    def get_simple_motif(motif):
        """
        Get the lexigraphically smallest repeat unit for a motif
        Arguments:
            motif: motif to get the repeat unit for
        Returns:
            lexigraphically smallest repeat unit for the motif
        """
        rus = SequenceAttributes.possible_ru(motif)
        return sorted(rus)[0]
    
    @staticmethod
    def get_terminal_label(row, upstream_max_dist=5, downstream_max_dist=5):
        # check if already labelled as external
        if row["internal_external"] == "external":
            if row["orientation"] == "+":
                if row["start"] < row["start_mei"]:
                    return "5'terminal"
                elif row["end"] > row["end_mei"]:
                    return "3'terminal"
            else: # reverse orientation -- flip what is 5' and 3'
                if row["end"] > row["end_mei"]:
                    return "5'terminal"
                elif row["start"] < row["start_mei"]:
                    return "3'terminal"
        # sort internal annotations by if they are within min_dist of the end or start of the MEI
        if row["orientation"] == "+":
            prime5_min_dist = min(abs(row["start"] - row["start_mei"]), abs(row["end"] - row["start_mei"]))
            prime3_min_dist = min(abs(row["start"] - row["end_mei"]), abs(row["end"] - row["end_mei"]))
            if prime5_min_dist <= upstream_max_dist:
                return "5'terminal"
            elif prime3_min_dist <= downstream_max_dist:
                return "3'terminal"
        elif row["orientation"] == "C":
            prime5_min_dist = min(abs(row["end"] - row["end_mei"]), abs(row["start"] - row["end_mei"]))
            prime3_min_dist = min(abs(row["end"] - row["start_mei"]), abs(row["start"] - row["start_mei"]))
            if prime5_min_dist <= upstream_max_dist:
                return "5'terminal"
            elif prime3_min_dist <= downstream_max_dist:
                return "3'terminal"
        return "internal"

class ComparisonTests: # TODO I could probaby make this less static and call from the df
    def __init__(self, mei_str_path=None, non_mei_str_path=None, additional_metadata_paths=None, duplicate_handling = "keep_all"):
        """
        The metadata fields are used to store which matadata has been added to the data or is to be added
        all metadata should be added to the data so it can be accessed by analysis functions
        Arguments:
            mei_str_path: path to the mei_str metadata
            non_mei_str_path: path to the non_mei_str metadata
            additional_metadata_paths: dictionary of additional metadata types and paths to load and preprocess
            duplicate_handling: how to handle duplicates in the data, default is to keep all, options: "keep_all", "drop_duplicated", "keep_first", "keep_last"
        
        """
        self.mei_str_path = mei_str_path
        self.non_mei_str_path = non_mei_str_path
        self.metadata = None
        self.tdb_data = None # data processed by TDBUtils
        self.longest_pure_segment_data = None # data loaded from LongestPureSegment processed data
        self.TRID_LocusID_Coords_map = None # stores a mapping of TRID to LocusID and the coordinates of each locus for going between data sources
        self.duplicate_handling = duplicate_handling # how to handle duplicates in the data, default is to keep all
        self.additional_metadata = additional_metadata_paths if additional_metadata_paths else defaultdict(list) # dictionary of additional metadata types and paths to load

        if mei_str_path and non_mei_str_path: # base needed input metadata
            self.metadata = self._load_metadata(mei_str_path, non_mei_str_path)
        # get locus mappings
        self.TRID_LocusID_Coords_map = self.get_locus_mappings()
    
    def describe_data(self, data_type = "all"):
        """
        Describe the data in the ComparisonTests object
        Arguments:
            data_type: dataframe(s) to describe, options: "all", "tdb_data", "longest_pure_segment_data"
        """
        if data_type == "all" or data_type == "tdb_data":
            print("TDB Data:")
            if self.tdb_data is not None:
                print(self.tdb_data.head())
                print(f"tdb_data has: {self.tdb_data.shape} dimensions")
            else:
                print("No TDB Data")
        if data_type == "all" or data_type == "longest_pure_segment_data":
            if self.longest_pure_segment_data is not None:
                print(f"Longest Pure Segment Data has {len(self.longest_pure_segment_data.keys())} methods")
                for agg_method in self.longest_pure_segment_data.keys():
                    print(f"Aggregation Method: {agg_method}")
                    print(self.longest_pure_segment_data[agg_method].head())
                    print(f"longest pure segment {agg_method} has: {self.longest_pure_segment_data[agg_method].shape} dimensions")
            else:
                print("No Longest Pure Segment Data")
        if self.metadata is not None:
            print("Metadata:")
            print(self.metadata.head())
            print(f"Metadata has: {self.metadata.shape} dimensions")
        if self.TRID_LocusID_Coords_map is not None:
            print("TRID to LocusID Mapping:")
            print(self.TRID_LocusID_Coords_map.head())
            print(f"TRID to LocusID Mapping has: {self.TRID_LocusID_Coords_map.shape} dimensions")
        if self.additional_metadata:
            print("Additional Metadata:")
            for metadata_type, metadata_paths in self.additional_metadata.items():
                print(f"additional metadata type: {metadata_type}:")
                print(metadata_paths)
        if self.duplicate_handling:
            print(f"Duplicate Handling: {self.duplicate_handling}")


    def get_locus_mappings(self):
        """
        Get the TRID to LocusID mappings from the tdb_data, longest_pure_segment_data, and metadata
        Sets self.TRID_LocusID_Coords_map to DataFrame with the TRID to LocusID to coordinate mappings
        """
        # TODO handle duplicate TRIDs in the mapping
        if self.metadata is None:
            raise ValueError("No metadata to get locus mappings from") 
        if self.TRID_LocusID_Coords_map is not None:
            return self.TRID_LocusID_Coords_map.drop_duplicates(subset = ["TRID","LocusID"])
        elif self.metadata is not None:
            # self.metadata has all the TRID to LocusID to coordinate mappings, this will just be a subset of that
            return self.metadata[["TRID", "LocusID", "chrom", "start", "end"]].drop_duplicates(subset = ["TRID","LocusID"])
        return None

    # TODO test this
    @staticmethod
    def _load_metadata(mei_str_path, non_mei_str_path):
        """
        Load the metadata for the mei_str and non_mei_str data
        Arguments:
            mei_str_path: path to the mei_str metadata
            non_mei_str_path: path to the non_mei_str metadata
        Returns:
            DataFrame: DataFrame with the metadata for the mei_str and non_mei_str data
        """
        mei_str_labelled = pd.read_csv(mei_str_path, sep="\t", header=None)
        non_mei_str_labelled = pd.read_csv(non_mei_str_path, sep="\t", header=None)

        mei_str_labelled.columns = ["chrom", "start_str", "end_str","str_id","str_motif","str_period","chrom_str",
                                    "start_mei","end_mei","mei_subfamily","num_merged","mei_family","dist","mei_str_id",
                                    "internal_external","bridge","start_labelled","end_labelled","TR_id","motif","motif_period","LocusID"]
        non_mei_str_labelled.columns = ["chrom", "start_str", "end_str","str_id","str_motif","str_period","chrom_str",
                                        "start_mei","end_mei","mei_subfamily","num_merged","mei_family","dist","chrom_labelled",
                                        "start_labelled","end_labelled","TR_id","motif","motif_period","LocusID"]
        
        mei_str_labelled['TRID'] = mei_str_labelled['str_id'].str.split('ID=', expand=True)[1].str.split(';', expand=True)[0]
        non_mei_str_labelled['TRID'] = non_mei_str_labelled['str_id'].str.split('ID=', expand=True)[1].str.split(';', expand=True)[0]

        mei_str_labelled["mei"] = True
        non_mei_str_labelled["mei"] = False

        all_metadata = pd.concat([mei_str_labelled[["LocusID","TRID", "chrom","start_str","end_str","str_period", "str_motif","mei_family","mei_subfamily","start_mei","end_mei","internal_external","bridge","mei"]], 
                                    non_mei_str_labelled[["LocusID","TRID", "chrom","start_str","end_str","str_period", "str_motif","mei"]]])
        # fill internal_external and mei_family with non_mei if they are not present
        all_metadata["internal_external"] = all_metadata["internal_external"].fillna("non_mei")
        all_metadata["mei_family"] = all_metadata["mei_family"].fillna("non_mei")
        all_metadata["mei_subfamily"] = all_metadata["mei_subfamily"].fillna("non_mei")
        all_metadata["bridge"] = all_metadata["bridge"].fillna("non_mei")

        # rename the start_str and end_str columns to match the other data
        all_metadata.rename(columns={"start_str":"start", "end_str":"end"}, inplace=True)
        return all_metadata
    
    def add_metadata(self, metadata_type, metadata_paths, data_type = "all", agg_methods = None):
        """
        Add additional metadata to the data DataFrame (either tdb_data or longest_pure_segment_data)
        All metadata should include the columns from the mapping of TRID to LocusID and coordinates 
        so they can be merged with the data for analysis regardless of source
        Arguments:
            metadata_type: type of metadata to add. Options: "satellite_data", "genomic_annotations"
            metadata_paths: dictionary of paths needed for the metadata to be loaded and preprocesssed
                            NOTE the keys of the dictionary should match the metadata_type
            data_type: dataframe(s) to add the metadata to, options: "all", "tdb_data", "longest_pure_segment_data"
            agg_methods: list of methods used to aggregate the longest_pure_segment_data, default is None
                        if data_type is "all" or "longest_pure_segment_data" and agg_method is None,
                        the metadata will be added to all agg_methods in longest_pure_segment_data dict

            NOTE I could probably be further refactored to take any data given the columns in the mapping it has 
        """
        if data_type not in ["all", "tdb_data", "longest_pure_segment_data"]:
            raise ValueError(f"Invalid data_type: {data_type}, options are: 'all', 'tdb_data', 'longest_pure_segment_data'")
        preprocessing_functions = {
            'satellite_data': self._preprocess_satellite_metadata,
            'genomic_annotations': self._preprocess_genomic_annotation_metadata,
            'exclusions': self._get_exclusions
            # Add more types and their corresponding preprocessing functions here
        }

        expected_path_keys = {
            'satellite_data': ["repeatmasker_hg38", "repeatmasker_all_alleles","allele_table"],
            'genomic_annotations': ["genomic_annotation"],
            'exclusions': [] # this may change, but this is the current possible inputs, really just needs to be a bed file
        }

        if metadata_type in preprocessing_functions:
            # check that the input metadata_paths keys correspond to the metadata_type preprocessing
            # expected keys for each metadata_type
            expected_keys = expected_path_keys[metadata_type]
            if set(metadata_paths.keys()) != set(expected_keys) and metadata_type != "exclusions":
                raise ValueError(f"Expected keys for metadata_type: {metadata_type} are: {expected_keys}, got: {metadata_paths.keys()}")
            additional_metadata = preprocessing_functions[metadata_type](metadata_paths)

            if data_type == "all" or data_type == "tdb_data":
                # check if the metadata has already been added
                if metadata_type == "satellite_data" and "satellite_hg38" in self.tdb_data.columns or metadata_type == "genomic_annotations" and  "genomic_annotation" in self.tdb_data.columns: # FIXME this needs to be generalized
                    print(f"Metadata: {metadata_type} already added to tdb_data!")
                    return #TODO figure out how I could have this skip to longest_pure_segment_data if it hasn't been added

                # Add the metadata to the tdb_data dataframe
                if self.tdb_data is not None:
                    orig_len = len(self.tdb_data)
                    self.tdb_data = pd.merge(self.tdb_data, additional_metadata, on=["LocusID","TRID","chrom","start","end"], how="left")
                    # assert the data has been added without changing the number of rows in data
                    try:
                        assert len(self.tdb_data) == orig_len
                        print(f"Metadata: {metadata_type} added successfully to tdb_data!", metadata_type)
                    except:
                        raise ValueError("Metadata did not merge correctly with tdb_data, before adding metadata: ", orig_len, " after adding metadata: ", len(self.tdb_data))
                else:
                    raise ValueError("No data to add metadata to")
                
            if data_type == "all" or data_type == "longest_pure_segment_data":
                if agg_methods is None:
                    agg_methods = list(self.longest_pure_segment_data.keys())
                # Add the metadata to the longest_pure_segment_data dataframe
                if self.longest_pure_segment_data is not None:
                    for agg_method in agg_methods:
                        if metadata_type == "satellite_data" and "satellite_hg38" in self.longest_pure_segment_data[agg_method].columns or metadata_type == "genomic_annotations" and  "genomic_annotation" in self.longest_pure_segment_data[agg_method].columns: # FIXME this needs to be generalized
                            print(f"Metadata: {metadata_type} already added to longest_pure_segment_data {agg_method}!")
                            continue # TODO see above about how we should be dealing with this

                        # TODO add case for genomic annotations after that is implemented
                        # drop the columns that came from the mapping that are not in longest_pure_segment_data
                        additional_metadata = additional_metadata.drop(columns=["LocusID","chrom","start","end"], errors="ignore")

                        orig_len = len(self.longest_pure_segment_data[agg_method])
                        self.longest_pure_segment_data[agg_method] = pd.merge(self.longest_pure_segment_data[agg_method], additional_metadata, on=["TRID"], how="left")
                        # assert the data has been added without changing the number of rows in data
                        try:
                            assert len(self.longest_pure_segment_data[agg_method]) == orig_len
                            print(f"Metadata: {metadata_type} added successfully to longest_pure_segment!", metadata_type)
                        except:
                            raise ValueError("Metadata did not merge correctly with longest_pure_segment_data, before adding metadata: ", orig_len, " after adding metadata: ", len(self.longest_pure_segment_data[agg_method]))
                else:
                    raise ValueError("No data to add metadata to")
    
    def _preprocess_satellite_metadata(self, metadata_paths):
        """
        Preprocess the satellite metadata
        Arguments:
            metadata_paths: dictionary of paths to the satellite metadata including repeatmasker 
                            satellites from hg38 and repeatmasker outputs from all alleles and allele_table
        Returns:
            DataFrame: DataFrame with the satellite metadata preprocessed

        NOTE this could be separated or at least the calculations could be helper functions
        """
        # Load the satellite metadata from the provided paths
        allele_table = pd.read_csv(metadata_paths["allele_table"])
        all_alleles_repeatmasker = pd.read_csv(metadata_paths["repeatmasker_all_alleles"], sep="\t", header=None)

        all_alleles_repeatmasker = all_alleles_repeatmasker.dropna()
        all_alleles_repeatmasker.columns = ["LocusID:allele","query_start","query_end","subfamily","family"]
        # split the LocusID:allele column
        all_alleles_repeatmasker["LocusID"] = all_alleles_repeatmasker["LocusID:allele"].str.split(":", expand=True)[0].astype(int)
        all_alleles_repeatmasker["allele_number"] = all_alleles_repeatmasker["LocusID:allele"].str.split(":", expand=True)[1].astype(int)
        allele_satellites = all_alleles_repeatmasker[all_alleles_repeatmasker["family"].str.contains("Satellite")]
        allele_satellites["satellite_allele"] = True

        # use the allele_table to get the number of total alleles per locus with a satellite
        allele_table = allele_table.merge(allele_satellites[["LocusID","allele_number","satellite_allele"]], on=["LocusID","allele_number"], how='left')
        allele_table["satellite_allele"] = allele_table["satellite_allele"].fillna(False)

        # calculate the percent of alleles labelled as satellite per locus 
        allele_table["percent_satellite_allele"] = allele_table.groupby("LocusID").satellite_allele.transform("sum")/allele_table.groupby("LocusID").allele_number.transform("count")
        # also add a column indicating at least one allele is labelled as satellite
        allele_table["satellite_allele_labelled"] = allele_table.groupby("LocusID").satellite_allele.transform("sum") > 0
        
        #TODO I changed the allele_table back to include all loci because when we merge back with the data we don't want nan values
        # reduce allele_table to unique on LocusID
        allele_table_locus = allele_table.drop_duplicates(subset = "LocusID")

        # load hg38 repeatmasker data
        hg38_repeatmasker = pd.read_csv(metadata_paths["repeatmasker_hg38"], sep="\t", header=None)
        hg38_repeatmasker.columns = ["chrom","start","end","repeat_subfamily","repeat_family"]
        hg38_repeatmasker["satellite_hg38"] = True

        # get the coordinates of our full input from the mapping table
        TRID_LocusID_Coords_map = self.TRID_LocusID_Coords_map.copy()
        # get the pyranges object for this data by renaming the chrom, start, and end columns to expected input
        hg38_repeatmasker_gr = pr.PyRanges(hg38_repeatmasker.rename(columns={"chrom":"Chromosome","start":"Start","end":"End"}))
        TRID_LocusID_Coords_map_gr = pr.PyRanges(TRID_LocusID_Coords_map.rename(columns={"chrom":"Chromosome","start":"Start","end":"End"}))
        overlaps = TRID_LocusID_Coords_map_gr.join(hg38_repeatmasker_gr, report_overlap= True, suffix="_repeatmasker")
        overlaps_df = overlaps.df
        overlaps_df["satellite_hg38"] = True
        overlaps_df = overlaps_df.rename(columns={"Chromosome":"chrom","Start":"start","End":"end"})
        # merge back to TRID_LocusID_Coords_map
        TRID_LocusID_Coords_map = pd.merge(TRID_LocusID_Coords_map, overlaps_df[["chrom","start","end","satellite_hg38"]], on=["chrom","start","end"], how="left")
        TRID_LocusID_Coords_map["satellite_hg38"] = TRID_LocusID_Coords_map["satellite_hg38"].fillna(False)

        # merge the allele_table locus with the TRID_LocusID_Coords_map on LocusID
        TRID_LocusID_Coords_map = pd.merge(TRID_LocusID_Coords_map, allele_table_locus[["LocusID","percent_satellite_allele","satellite_allele_labelled"]], on="LocusID", how="left")
        TRID_LocusID_Coords_map["percent_satellite_allele"] = TRID_LocusID_Coords_map["percent_satellite_allele"].fillna(0)
        TRID_LocusID_Coords_map["satellite_allele_labelled"] = TRID_LocusID_Coords_map["satellite_allele_labelled"].fillna(False)
        
        return TRID_LocusID_Coords_map.drop_duplicates(subset = "LocusID")

    def _get_exclusions(self,exclusion_paths):
        """
        Get the exclusions from the exclusion paths
        Arguments:
            exclusion_paths: dictionary of paths to the exclusion files (ie centromeres, telomeres, blacklist, satellites ect.)
        Returns:
            DataFrame: DataFrame with the exclusions
        """
        exclusions = []
        for exclusion_type, exclusion_path in exclusion_paths.items():
            exclusion = pd.read_csv(exclusion_path, sep="\t", header=None)[[0,1,2]] # only get coordinate columns
            exclusion.columns = ["Chromosome","Start","End"]
            exclusion["exclusion_type"] = exclusion_type
            exclusions.append(exclusion)
        exclusions_df = pd.concat(exclusions)
        # merge with our mapping table
        TRID_LocusID_Coords_map = self.TRID_LocusID_Coords_map.copy()
        TRID_LocusID_Coords_map_gr = pr.PyRanges(TRID_LocusID_Coords_map.rename(columns={"chrom":"Chromosome","start":"Start","end":"End"}))
        exclusions_gr = pr.PyRanges(exclusions_df)
        overlaps = TRID_LocusID_Coords_map_gr.join(exclusions_gr, report_overlap= True, suffix="_exclusion")
        overlaps_df = overlaps.df
        # merge back to TRID_LocusID_Coords_map
        overlaps_df = overlaps_df.rename(columns={"Chromosome":"chrom","Start":"start","End":"end"})
        TRID_LocusID_Coords_map = pd.merge(TRID_LocusID_Coords_map, overlaps_df[["chrom","start","end","exclusion_type"]], on=["chrom","start","end"], how="left")
        TRID_LocusID_Coords_map["exclusion_type"] = TRID_LocusID_Coords_map["exclusion_type"].fillna("keep")
        return TRID_LocusID_Coords_map.drop_duplicates(subset = "LocusID")

    def _preprocess_genomic_annotation_metadata(self,metadata_path_dict):
        '''
        Preprocess the genomic annotation metadata
        Arguments:
            metadata_path_dict: a dictionary with key set to genomic_annotation and value path to the genomic annotation bed file
        '''
        metadata_path = metadata_path_dict["genomic_annotation"]
        metadata= pd.read_csv(metadata_path, sep="\t", header=None)

        metadata.columns = ["Chromosome","Start","End","genomic_annotation"]

        #Lets try to merge all the metadata into one pyranges object then intersect with the TRID mapping dataframe
        metadata_gr = pr.PyRanges(metadata)
        TRID_LocusID_Coords_map = self.TRID_LocusID_Coords_map.copy()
        TRID_LocusID_Coords_map = TRID_LocusID_Coords_map.rename(columns={"chrom":"Chromosome","start":"Start","end":"End"})
        TRID_LocusID_Coords_map_gr = pr.PyRanges(TRID_LocusID_Coords_map[["Chromosome","Start","End"]])

        # get the overlaps
        overlaps = TRID_LocusID_Coords_map_gr.join(metadata_gr, report_overlap=True, suffix="_metadata") # there are duplicates in the TRID map coordinates because of multiple motif annotations
        overlaps_df = overlaps.df.sort_values(by = "Overlap",ascending=False).drop_duplicates(subset=["Chromosome","Start","End"],keep="first") # only keep one annotation per locus in our TR set

        # merge back to the TRID_LocusID_Coords_map
        TRID_LocusID_Coords_map = TRID_LocusID_Coords_map.merge(overlaps_df[["Chromosome","Start","End","genomic_annotation"]], on=["Chromosome","Start","End"], how='left')
        # rename the columns back to the original
        TRID_LocusID_Coords_map = TRID_LocusID_Coords_map.rename(columns={"Chromosome":"chrom","Start":"start","End":"end"})

        return TRID_LocusID_Coords_map.drop_duplicates(subset = "LocusID")


    def _add_mei_labels(self, data):
        """
        Add mei labels to the data
        Arguments:
            data: DataFrame to add mei labels to
        Returns:
            DataFrame: DataFrame with mei labels added
        """
        data["Alu"] = data["mei_family"].str.contains("Alu")
        data["L1"] = data["mei_family"].str.contains("L1")
        data["SVA"] = data["mei_family"].str.contains("SVA")
        data["Alu_subfamily"] = data["mei_subfamily"].apply(lambda x: "AluY" if "AluY" in x else ("AluS" if "AluS" in x else ("AluJ" if "AluJ" in x else "non_Alu")))
        data.loc[~data["mei"], "Alu_subfamily"] = "non_mei"
        return data
    
    def _duplicate_handling(self, data, duplicate_column):
        """
        Handle duplicates in the data based on the duplicate_handling attribute
        Arguments:
            data: DataFrame to handle duplicates in
            duplicate_column: column to check for duplicates
        Returns:
            DataFrame: DataFrame with duplicates handled
        """
        if self.duplicate_handling == "drop_duplicated":
            return data[~data.duplicated(subset = duplicate_column, keep=False)] # drops all loci that appear more than once
        elif self.duplicate_handling == "keep_first":
            return data.drop_duplicates(subset = duplicate_column,keep="first")
        elif self.duplicate_handling == "keep_last":
            return data.drop_duplicates(subset = duplicate_column,keep="last")
        
        return data

    #TODO test this since updated to use the metadata function
    def load_length_attributes(self, cleaned_data = None, observed_length_path = None, 
                                length_polymorphism_path = None):
        """
        Load the length attributes for all data and add add labels to the data for downstream analysis.
        If the length attributes have already beed calculated, load them from the provided paths,
        otherwise calculate them and save them to the provided paths.

        NOTE I have moved the initial metadata path inputs to constructor and if additional metadata is included 
        on initialization, the dictionary needed to load it is constructed and add_metadata can be called here
        to initialize the data with all metadata if desired

        Arguments:
            cleaned_data: tdb database object to conduct the analysis on, 
                           None if the length attributes have already been calculated
            observed_length_path: path to the saved observed length data if cleaned_data is None
                                   path to save the observed length data if cleaned_data is not None            
            length_polymorphism_path: path to the length polymorphism data if cleaned_data is None
                                      path to save the length polymorphism data if cleaned_data is not None
            LongestPureSegment_path: path to the LongestPureSegment data if it exists.
            LongestPureSegment_methods: list of methods to use to calculate the longest pure segment per locus
                                        options: "greatest_N_motifs", "weighted_average_LPS", "naiive_drop_duplicates"
                                        default: ["greatest_N_motifs"]

        Sets the tdb_data attribute to the data with the length attributes added

        """
        if cleaned_data is None and os.path.exists(observed_length_path) and os.path.exists(length_polymorphism_path):
            # the length attributes have already been calculated, load them
            print("Loading precalculated length data...")
            observed_length = pd.read_csv(observed_length_path)
            length_polymorphism = pd.read_csv(length_polymorphism_path)
            print("Length attributes loaded successfully!")
        else:
            # Calculate the length attributes
            try:
                observed_length = TDBUtils.get_observed_allele_lengths(cleaned_data, cleaned_data['sample'].keys())
                length_polymorphism = TDBUtils.length_polymorphism_score(cleaned_data, min_af = 0.0)

                # Save the length attributes
                observed_length.to_csv(observed_length_path, index=True) # the above functions change the index to LocusID
                length_polymorphism.to_csv(length_polymorphism_path,index=True)
            except:
                print("Error calculating length attributes from given cleaned_data, tdb_data will remain null")
        
        print("Length attributes loaded successfully! Calculating stats...")
        # calculate summary stats from observed_length
        observed_lengths_stats = observed_length.groupby("LocusID")["allele_length"].agg(["mean", "median", "std"]).reset_index()

        print("Stats calculated successfully! Combining data...")
        #Combine length metrics calculalated from the cleaned data
        all_length_data = length_polymorphism.merge(observed_lengths_stats, on="LocusID", how="inner")

        # add initial metadata to the data
        all_length_data = pd.merge(left = all_length_data, right = self.metadata, how = "left", on = ["LocusID","chrom","start","end"]) # FIXME I think this is all I need
        print("Metadata added successfully!")

        # add Alu, L1, and SVA labels
        self.tdb_data = self._add_mei_labels(all_length_data) #NOTE this MUST be set before adding additional metadata, cannot be null
        print("MEI family labels added successfully!")

        # add any initial additional metadata
        for metadata_type, metadata_paths in self.additional_metadata.items(): # should not add anything if there is no additional metadata
            self.add_metadata(metadata_type, metadata_paths, data_type="tdb_data")
        
        # handel duplicates as designated in the constructor
        print(f"Handling duplicates with {self.duplicate_handling}, original dimensions: ", str(self.tdb_data.shape[0]))
        self.tdb_data = self._duplicate_handling(self.tdb_data, "LocusID")
        print("Duplicates handled successfully! New dimensions: ", str(self.tdb_data.shape[0]))


    # TODO test this since refactored
    def load_longest_pure_segment(self,LongestPureSegment_path, LongestPureSegment_methods = ["greatest_N_motifs"]):
        """
        Load the LongestPureSegment data and aggregate it by the methods defined in the function

        Arguments:
            LongestPureSegment_path: path to the LongestPureSegment data
            LongestPureSegment_methods: list of methods to use to calculate the longest pure segment per locus
                                        options: "greatest_N_motifs", "weighted_average_LPS", "naiive_drop_duplicates"
                                        default: ["greatest_N_motifs"]
        set longest_pure_segment_data: DataFrame with the LongestPureSegment data aggregated by the methods defined in the function
        """
        try:
            LongestPureSegment = pd.read_csv(LongestPureSegment_path, sep="\t")
            
            # aggregate the LongestPureSegment data by the methods defined in the function
            LPS_aggregated = ComparisonTests._LongestPureSegment_aggregate(LongestPureSegment, LongestPureSegment_methods)
            print("LongestPureSegment aggregated successfully! Adding metadata for groups: ", LPS_aggregated.keys())

            print("LongestPureSegment aggregated successfully! Adding metadata...")
            # initialize longest_pure_segment as a dictionary of DataFrames
            self.longest_pure_segment_data = {}

            for agg_method, LPS_data in LPS_aggregated.items():
                # add mei labels
                print("Adding labels to: ", agg_method)
                # print(LPS_data.head())
                # subset to only the TRIDs that are in the mei_str_labelled or non_mei_str_labelled
                LPS_data = LPS_data[LPS_data["TRID"].isin(self.metadata["TRID"].tolist())]
                LPS_data['mei'] = LPS_data['TRID'].isin(self.metadata[self.metadata.mei]['TRID'].tolist())
                print("MEI labels added successfully!")
                # add internal external labels
                LPS_data = pd.merge(left = LPS_data, right = self.metadata[["TRID", "internal_external"]], how = "left", on = "TRID")
                LPS_data["internal_external"] = LPS_data["internal_external"].fillna("non_mei") # FIXME probably not necessary anymore
                print("Internal external labels added successfully!")

                # TODO figure out how we want to add period and motif to this data given that it is calculated separately and doesn't aggregate well
                LPS_data = pd.merge(left = LPS_data, right = self.metadata[["TRID","mei_family","mei_subfamily","str_motif","str_period"]], how = "left", on = "TRID")
                print("MEI family labels added successfully!")

                self.longest_pure_segment_data[agg_method] = self._add_mei_labels(LPS_data) #NOTE this MUST be set before adding additional metadata, cannot be null
                print("MEI family labels added successfully!")

                # add any initial additional metadata
                for metadata_type, metadata_paths in self.additional_metadata.items(): # should not add anything if there is no additional metadata
                    self.add_metadata(metadata_type, metadata_paths, data_type="longest_pure_segment_data", agg_method = agg_method)
                    print(f"Additional metadata {metadata_type} added successfully!")

                # add the Alu subfamily labels -- FIXME this is not working, try to do this after we've output the data
                #LPS_data["Alu_subfamily"] = LPS_data["mei_subfamily"].apply(lambda x: "AluY" if "AluY" in x else ("AluS" if "AluS" in x else ("AluJ" if "AluJ" in x else "non_Alu")))
                #LPS_data.loc[~LPS_data["mei"], "Alu_subfamily"] = "non_mei"

                LPS_aggregated[agg_method] = LPS_data

        except:
            print("WARNING: LengthPolymorphismScore does not exist or something went wrong, returning None")


    @staticmethod
    def _weighted_avg_and_std(df, weight_col, mean_col, std_col):
        """
        Calculate the weighted average and standard deviation.

        Parameters:
        df (DataFrame): DataFrame containing the data.
        weight_col (str): Column name for the weights.
        mean_col (str): Column name for the means.
        std_col (str): Column name for the standard deviations.

        Returns:
        DataFrame: DataFrame with weighted average and standard deviation for each TRID.
        """
        # Calculate the weighted mean
        df['weighted_mean'] = df[mean_col] * df[weight_col]
        weighted_mean = df.groupby('TRID')['weighted_mean'].sum() / df.groupby('TRID')[weight_col].sum()
        weighted_mean.name = "weighted_mean"
        df.drop(columns=["weighted_mean"], inplace=True)


        # Broadcast the weighted mean to match the dimensions of the DataFrame
        df = df.merge(weighted_mean.reset_index(), on='TRID')

        # Calculate the weighted variance
        df['weighted_variance'] = df[weight_col] * (df[std_col]**2 + (df[mean_col] - df['weighted_mean'])**2)
        weighted_variance = df.groupby('TRID')['weighted_variance'].sum() / df.groupby('TRID')[weight_col].sum()

        # Calculate the weighted standard deviation
        weighted_std = np.sqrt(weighted_variance)

        # Combine the results into a DataFrame
        result = pd.DataFrame({
            'Mean': weighted_mean,
            'Stdev': weighted_std
        }).reset_index()

        return result

    @staticmethod # NOTE when we are doing thing for length analysis we will not always maintain N_motifs
    def _LongestPureSegment_aggregate(data,methods):
        """
        Calculate the longest pure segment for each locus by an established method
        Arguments:
            data: LongestPureSegment data already loaded from tsv
            methods: list of methods to use to calculate the longest pure segment per locus
                     options: "greatest_N_motifs", "weighted_average_LPS", "naiive_drop_duplicates"
        Returns:
            dictionary of DataFrame(s) with TRID and LongestPureSegment data with each row being a unique locus
        """
        LPS_aggregated = {}
        #drop na values
        data = data.dropna()
        # method 1: sort by N_motifs and keep the row with the highest N_motifs as the most representative
        if "greatest_N_motifs" in methods:
            data_greatest_N_motifs = data.sort_values("N_motif", ascending=False).drop_duplicates("TRID")
            LPS_aggregated["greatest_N_motifs"] = data_greatest_N_motifs
        if "weighted_average_LPS" in methods:
            data_weighted_avg_LPS = ComparisonTests._weighted_avg_and_std(data, "N_motif", "Mean", "Stdev")
            LPS_aggregated["weighted_average_LPS"] = data_weighted_avg_LPS
        if "naiive_drop_duplicates" in methods:
            data_naive_drop_duplicates = data.drop_duplicates("TRID")
            LPS_aggregated["naive_drop_duplicates"] = data_naive_drop_duplicates
        if not methods: # keep all including duplicates
            LPS_aggregated["all"] = data
        return LPS_aggregated

    @staticmethod
    def add_satellite_labels(hg38_satellite_path, allele_repeatmasker):
        # TODO get this from the notebook
        return
    @staticmethod # TODO add compatibility for saving to a file instead of plotting to console
    def plot_histogram_comparisons(groups, data, labels = {},log=True, bins=100): # This will be loaded with tdb utils in the future but I dont want to restart the kernel
        """
        Plot histograms of the data for each group
        Arguments:
            groups: list of group names
            data: dictionary of dataframes with group names as keys
                    subset of dataframes already done when passed to this function
            labels: dictionary of labels for the groups including the metric per group
            log: boolean, whether to plot the data on a yscale = log
        """
        for group in groups:
            plt.hist(data[group], bins=bins)
            if log:
                plt.yscale('log')
            plt.axvline(data[group].mean(), color='r', label='mean: '+str(data[group].mean()))
            plt.axvline(data[group].median(), color='g', label='median' + str(data[group].median()))
            plt.legend()
            if not labels:
                plt.title('Histogram of ' + group + ", n = " + str(len(data[group])))
                plt.xlabel(group)
            else:
                plt.title('Histogram of ' + labels[group] + ", n = " + str(len(data[group])))
                plt.xlabel(labels[group])
            plt.ylabel('Frequency')
            plt.show()
            # plt.gca().clear()

        return

    @staticmethod
    def manwhitneyu_test(group1, group2, metric_col, labels = {"group1":"group1","group2":"group2"}): 
        """
        Perform a Mann-Whitney U test on two groups of data
        Arguments:
            group1: data from first group to compare
            group2: data from second group to compare
            metric_col: column name for the metric to compare

            labels: dictionary of labels for the groups
        Returns:
            p-value from the Mann-Whitney U test for the two groups of variables
        """
        print(f'Comparing {metric_col} across groups {labels["group1"]} and {labels["group2"]}')
        stat, all_pval = mannwhitneyu(group1[metric_col], group2[metric_col])
        print('All data: n = ', len(group1) + len(group2))
        print(f'{labels["group1"]} data: n1 = ', len(group1))
        print(f'{labels["group2"]} data: n2 = ', len(group2))
        print('Statistics=%.3f, p=%.3f' % (stat, all_pval)) # p-value is significant -- distributions are different but are the means of this distribution a good metric to compare?
        print()

        return all_pval

    @staticmethod
    def kruskalwallis_test(groups, metric_col, variable_groups=None,labels=None):
        """
        Perform a Kruskal-Wallis H test on multiple groups of data
        Arguments:
            groups: list of dataframes for each group
            metric_col: column name for the metric to compare
            labels: list in the same order as groups of labels for the groups
        Returns:
            p-value from the Kruskal-Wallis H test for the groups of variables
        """
        grouped_metric = [group[metric_col] for group in groups]
        stat, p = kruskal(*grouped_metric)
        
        print('All data: n = ', len(pd.concat(groups)))
        for i, group in enumerate(groups):
            if labels is not None:
                print(f'{labels[i]} data: n = ', len(group))
            else:
                print(f'Group {i} data: n = ', len(group))
        print('Statistics=%.3f, p=%.3f' % (stat, p))
        print()
        
        return p

    @staticmethod
    def get_percent_variable(data, variable_col, invariable_value, group_col):
        """
        Get the percentage of a variable that is a certain value
        Arguments:
            variable_col: column of the variable to analyze
            invariable_value: value to check for the percentage of
            group_col: column of the group to analyze
        Returns:
            percentage of the variable that is the invariable_value
        """
        return data[data[variable_col] == invariable_value].groupby(group_col).size() / data.groupby(group_col).size() * 100

    @staticmethod
    def plot_violin_comparison(groups, data, labels = {}): #FIXME need to figure out how we want to make this generalizable
        """
        Plot violin plots of the data for each group
        Arguments:
            groups: list of group names
            data: dictionary of dataframes with group names as keys
                  subset of dataframes already done when passed to this function
            labels: dictionary of labels for the groups
        """
        for group in groups:
            sb.violinplot(data=data[group], inner='point')
            plt.title('Violin plot of ' + labels.get(group, group))
            plt.show()
            plt.gca().clear()

        return
    @staticmethod
    def plot_boxplot_multi_comparison(hue_col_name, data, x_col, y_col, labels = {}, showfliers=True, log_scale=False, horizontal=True, flierprops = dict(marker='o', markeredgecolor = 'grey',alpha=0.25, markersize=5)):
        """
        Plot boxplots to compare the quantiles of the data for each group using seaborn -- assumes multiple groups
        Arguments:
            hue_col_name: column name for the groups
            data: DataFrame with the data
            x_col: column name for the x-axis
            y_col: column name for the y-axis
            labels: dictionary of labels for the plot including title, xlabel, ylabel, and legend_title
            showfliers: boolean, whether to show the outliers
            log_scale: boolean, whether to plot the data on a log scale
            horizontal: boolean, whether to plot the boxplots horizontally
            flierprops: dictionary of properties for the fliers

        """
        plt.figure(figsize=(12, 8))
        sb.boxplot(x=x_col, y=y_col, hue=hue_col_name, data=data, showfliers=showfliers, flierprops=flierprops)
        plt.title(labels["title"])
        plt.xlabel(labels["xlabel"])
        plt.ylabel(labels["ylabel"])
        if log_scale:
            plt.yscale('log')
        plt.legend(title= labels["legend_title"])
        plt.show()

        # Display the number of rows per group and print the median, mean, and standard deviation of each group
        grouped = data.groupby([x_col,hue_col_name])[y_col]
        for name, group in grouped:
            print(f"Group: {str(name[0])}, {str(name[1])}")
            print(f"Number of rows: {len(group)}")
            print(f"Median: {group.median()}")
            print(f"Mean: {group.mean()}")
            print(f"Standard Deviation: {group.std()}")
            print()
        return
    @staticmethod
    def plot_boxplot_comparison(groups, data, labels = {}, showfliers=True, log_scale=True, horizontal=True, flierprops = dict(marker='o', markeredgecolor = 'grey',alpha=0.25, markersize=5)):
        """
        Plot boxplots to compare the quantiles of the data for each group
        Arguments:
            groups: list of group names
            data: dictionary of dataframes with group names as keys
                    subset of dataframes already done when passed to this function
            labels: dictionary of labels for the groups and title, title should be "title" key
            showfliers: boolean, whether to show outliers
            log_scale: boolean, whether to plot the data on a log scale
            horizontal: boolean, whether to plot the boxplot horizontally
            flierprops: dictionary of properties for the fliers (outlier markers)
        """
        data_to_plot = [data[group] for group in groups]
        group_sizes = [len(data[group]) for group in groups]
        labels_with_size = [f"{labels.get(group, group)}\nn = {size}" for group, size in zip(groups, group_sizes)]
        plt.boxplot(data_to_plot, labels=[labels.get(group, group) for group in groups], showfliers=showfliers, vert=not horizontal, flierprops=flierprops)
        if "title" in labels:
            plt.title(labels["title"])
        else:
            plt.title('Comparison of ' + ', '.join(groups))
        if horizontal:
            plt.xlabel('Values')
            plt.ylabel('Groups')
            plt.yticks(ticks=range(1, len(groups) + 1), labels=labels_with_size)

        else:
            plt.xlabel('Groups')
            plt.ylabel('Values')
            plt.xticks(ticks=range(1, len(groups) + 1), labels=labels_with_size)
        if log_scale:
            if horizontal:
                plt.xscale('log')
            else:
                plt.yscale('log')
        plt.show()

        for group in data.keys():
            print(group, "mean:", data[group].mean(), "median:", data[group].median(), "std:", data[group].std())

        return

    @staticmethod
    def pairwise_histogram(data, x_col,hue_col, group_col, labels,log = True):
        """
        Plot histograms of the data for each group on the 
        Arguments:
            data: DataFrame with the data
            x_col: column name for the x-axis metric to plot
            hue_col: column name for the groups to be plotted on the same histogram
            group_col: column name for the groups to be plotted separately
            labels: dictionary of labels for the plot including title, xlabel, ylabel
            log: boolean, whether to plot the data on a log scale
        """
        for group in data[group_col].unique():
            hues = data[data[group_col] == group][hue_col].unique()
            # get range of data in this group to set bins for the histogram
            min_val = data[data[group_col] == group][x_col].min()
            max_val = data[data[group_col] == group][x_col].max()
            bins = np.linspace(min_val, max_val, 100)
            for hue in hues:
                plt.hist(data[data[group_col] == group][data[hue_col] == hue][x_col], label = hue, alpha = 0.5, bins=bins)
            if log:
                plt.yscale('log')
            plt.title(labels["title"] + ": " + str(group))
            plt.xlabel(labels["xlabel"])
            plt.ylabel(labels["ylabel"])
            plt.legend(title = hue_col)
            plt.show()

        return

    @staticmethod
    def variability_barplot(variability_dict, var_name, hue,title):
        # Transform the dictionary to a dataframe
        df = pd.DataFrame(variability_dict).reset_index().melt(id_vars=hue, var_name=var_name, value_name='Percent Variable')

        # Plot the boxplot
        plt.figure(figsize=(12, 6))
        sb.barplot(x=var_name, y='Percent Variable', hue=hue, data=df)
        plt.title(title)
        plt.xlabel(var_name)
        plt.ylabel('Percent Variable')
        plt.legend(title=hue)
        plt.show()

class DivergenceUtils():
    '''
    Helper functions for processing and plotting divergence
    and related data
    '''
    @staticmethod
    def divergence_support_filter(data, min_regions = 50):
        '''
        given a dataframe with divergence already included
        filter rows with divergence measures that are supported by
        less than min_regions
        
        Arguments:
            data: DataFrame with the data
            min_regions: minimum number of regions to support the divergence measure
        Returns: 
            DataFrame: DataFrame with the filtered data
        '''
        data_counts = data['divergence'].value_counts()
        data_filtered = data[data['divergence'].isin(data_counts[data_counts >= min_regions].index)]
        return data_filtered
    @staticmethod
    def plot_alu_family_dist(data):
        '''
        given the data with different Alu families, plot the frequency of
        each family as a function of divergence
        Arguments:
            data: DataFrame with the data including mei_family, divergence, and mei_subfamily 
                  columns
        Returns:
                nothing, plots the distribution of Alu families
        '''
        data = data[data['mei_family'] == 'SINE_Alu']

        # Group mei_subfamilies into AluY, AluS, and AluJ
        data['group'] = data['mei_subfamily'].apply(lambda x: 'AluY' if 'AluY' in x else ('AluS' if 'AluS' in x else ('AluJ' if 'AluJ' in x else 'Other')))

        # Get unique groups
        alu_groups = data['group'].unique()
        plt.figure(figsize=(12, 6))
        for group in alu_groups:
            if group == 'Other':
                continue
            group_mei_str = data[data['group'] == group]
            sb.histplot(data=group_mei_str, x='divergence', label=group, alpha=0.5, bins=100)
            plt.title('Distribution of Divergence for Alu Groups for Alus with STRs')
            plt.xlabel('Divergence')
            plt.ylabel('Count')
            plt.legend()
            plt.show()
        return
    
    @staticmethod
    def plot_rolling_proportion(data, groups, window_size = 15):
        '''
        given a column with groups to plot on the same plot, plot the metric
        as a function of the rolling mean of the metric as a function of divergence
        Arguments:
            data: DataFrame with the data including divergence, and the columns with the groups and metric
            groups: list of groups to plot
            window_size: size of the rolling window
        Returns:
            to_plot: the grouped dataframe that was plot in case this is wanted later
        '''
        data_proportions = data.groupby(['divergence', groups]).size().unstack().fillna(0)
        data_proportions = data_proportions.div(data_proportions.sum(axis=1), axis=0).rolling(window=window_size).mean()
        plt.figure(figsize=(12, 6))
        for var in data_proportions.columns:
            sb.lineplot(data=data_proportions, x=data_proportions.index, y=var, label=f'{var}')
        plt.title(f'Rolling Proportion of {groups} by Divergence, window size = {window_size}')
        plt.xlabel('Divergence')
        plt.ylabel(f'Proportion of {groups}')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.show()
        return data_proportions
    
    @staticmethod
    def plot_metric_by_divergence(data, groups, metric, window_size = 15):
        '''
        Plots a given metric in data as a function of divergence for groups in groups column
        '''
        # get the mean of the metric per group per divergence
        data_means_grouped = data.groupby(['divergence',groups])[metric].mean().reset_index()
        # apply smoothing as the rolling mean for windows of window_size
        data_means_grouped[metric] = data_means_grouped[metric].rolling(window=window_size).mean()
        plt.figure(figsize=(12, 6))
        for group in data_means_grouped[groups].unique():
            sb.lineplot(data=data_means_grouped[data_means_grouped[groups] == group], x='divergence', y=metric, label=f'{group}')
        plt.title(f'{metric} by Divergence, window size = {window_size}')
        plt.xlabel('Divergence')
        plt.ylabel(f'{metric}')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.show()
        return data_means_grouped



    
