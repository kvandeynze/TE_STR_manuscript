#!/bin/bash
# we need to figure out (1) how to deal with commas in the motifs for these (2) how the heck does TRGT deal with reporting lengths for multiple motifs --
# the inputs place a large window around some of these such that it really wouldn't be fair to compare the lengths of regions defined this way if they
# represent multiple neighboring TRs that are all polymorphic

# daTR regions that have been filtered out because of these issues
# HD and DM2 -- these were defined as a VC and the TRGT input includes a large window around the TRs to include even non overlapping TRs 
# SCA27B and FXS -- These include interruption sequence in their definition which I don't disagree with because they are within the same bounds
# ie the definition of the region is the same as the definition of the daTRs but it include subregions where there may be another sequence -- fine
# There are also multiple definitions for some of the same regions where the periods are different or the motifs are different -- these got filtered out 
# when I check for overlaps but if we want to still keep them we need to figure out a fool proof way to be consistently choosing the definition we want
# see a few of the FAME loci where one is defined as both a tetranucleotide and a pentanucleotide repeat 
# and the other includes a homopolymer definition upstream that overlaps the pentanucleotide definition

# Variation Cluster regions as defined by the authors -- A lot of the regions that overlap with the daTRs are not well matched and don't have the expected motifs
# some look separaple like SCA8 and HD
# some of these overlap but aren't the same STR (DM1, BSS, ...) presumably the acutal daTR is defined in the isolated set not the VC and was dropped either when I dropped overlaps
# with VCs or when I dropped overlaping regions
VC_regions=/data/projects/nanopore/RepeatExpansion/coordinates/variation_clusters_v1.hg38.TRGT.bed
ALL_regions=/data/projects/nanopore/RepeatExpansion/coordinates/variation_clusters_and_isolated_TRs_v1.hg38.TRGT.bed

# some of these are defined in VCs and isolated TRs and others are only defined as one or the other --
# if there is an isolated definition we want to take that one over the VC deifnition because it is easier to parse and more specific
# for instance, if you intersect VCs and VCs and isolated TRs there are nearly twice as many regions as in the isolated TRs corresponding to the VC ranges
# There are 28 daTRs that intersect VCs, of these 19 of them have multiple entries in the varaition cluster and isolated file, hence these would be discarded
# in the isolated file and would also be discarded if we discardoverlaps since these are duplicate annotations

# output files
OUT_DIR=/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/sorted_VCs_overlaps
VCs_merged=${OUT_DIR}/variation_clusters_merged.bed
VC_independent=${OUT_DIR}/variation_clusters_independent.bed # regions that are independent of other VCs
VC_only=${OUT_DIR}/variation_clusters_only.bed # regions that are only found in varaition clusters
ISOLATED_only=${OUT_DIR}/isolated_TRs_only.bed # regions that are only found in isolated TRs that do not overlap any other annotations, can now add back bookend to this
BOTH=${OUT_DIR}/both.bed # regions that have annotations in both VCs and isolated TRs
BOTH_PAIRS=${OUT_DIR}/both_pairs.bed # regions that have annotations in both VCs and isolated TRs and are overlapping
BOOKEND_VC=${OUT_DIR}/bookend_vc.bed # regions that are bookend to other annotations in VCs
OVERLAP_no_VC_pairs=${OUT_DIR}/overlap_no_vc_pairs.bed # regions that overlap with other annotations in isolated TRs
BOOKEND_no_VC=${OUT_DIR}/bookend_no_vc.bed # regions that are bookend to other annotations in isolated TRs
OVERLAP_no_VC=${OUT_DIR}/overlap_no_vc.bed 

# get the regions that are only in the VCs
# VC_regions is a subset of ALL_regions -- this set is the regions that are unique to the VC set
# this means they will only come up once in the ALL_regions file

#merge the VC regions to get unique regions with respect to VCs (there are ~100 overlapping VCs)
bedtools merge -c 4 -o collapse,count -i $VC_regions > $VCs_merged
# get independent VCs 
awk '$5 == 1' $VCs_merged | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > $VC_independent # drops the 11 regions with multiple annotations in VC file
# get the regions that are only in the VCs
bedtools intersect -b $VC_independent -a $ALL_regions -wa| bedtools sort -i - > tmp # -wa should save the original entry in ALL_regions which has both VC and non-VC annotations -- ensures bedtools doesn't mess with the coordinates
# merge 
bedtools merge -c 4 -o collapse,count -i tmp > tmp2
# drop all regions that have count more than 1 by intersecting with the original VC file
awk '$5 == 1' tmp2 | bedtools intersect -a $VC_independent -b - -wa > $VC_only # about half of VCs are unique to the VC set (142K unique vs 130k with multiple annotations)

# get the regions that are defined in both VCs and isolated TRs -- take all definitions for these so we can figure out what to do with them
# To Note: there are regions here that have overlapping annotations in VC and isolated according to the coordinates but they aren't overlapping according to the TRIDs or overlap very little -- may be worth quatifying this
awk '$5 > 1' tmp2 | bedtools intersect -a $ALL_regions -b - -wa > $BOTH # since there are overlapping/duplicate annotations in this set, there are 281K entries
# move the merged regions to the overlap VC overlap file so we can look at pairs of overlapping regions
awk '$5 > 1 {print}' tmp2 > $BOTH_PAIRS # there are 130K pairs of overlapping regions that are in both VC and isolated TRs

# get bookend VCs and isolated from this set to separate
bedtools intersect -C -a $BOTH -b $BOTH | awk '$NF == 1 {print}' > $BOOKEND_VC # ~10% if the overlapping regions
# take these out of the overlap_no_vc file
bedtools intersect -v -a $BOTH -b $BOOKEND_VC > tmp5
mv tmp5 $BOTH

# purely isolated TRs
bedtools intersect  -v -a $ALL_regions -b $VC_only | bedtools intersect -v -a - -b $BOTH | bedtools intersect -v -a - -b $VCs_merged > tmp3 # 4116348 regions

# divide into isolated only and overlapping regions
bedtools merge -c 4 -o collapse,count -i tmp3 > tmp4
awk '$5 == 1' tmp4 | bedtools intersect -a tmp3 -b - -wa > $ISOLATED_only # 3.8M regions, this is what we had previously so this is good

awk '$5 > 1' tmp4 | bedtools intersect -a tmp3 -b - -wa > $OVERLAP_no_VC # 312k regions have overlaps, this is what we had previously so this is good
# get all regions in overlaps that are bookend and not overlapping
bedtools intersect -C -a $OVERLAP_no_VC -b $OVERLAP_no_VC | awk '$NF == 1 {print}' > $BOOKEND_no_VC # ~10% if the overlapping regions
# take these out of the overlap_no_vc file
bedtools intersect -v -a $OVERLAP_no_VC -b $BOOKEND_no_VC > tmp5
mv tmp5 $OVERLAP_no_VC

awk '$5 > 1 {print}' tmp4 > $OVERLAP_no_VC_pairs 

# FINAL NUMBERS             TOTAL                     daTRs in set
# VC_only:                  142,752                   8: SCA8, SCA27B, SCA31, RCPS, SCA36, DM2, HD, FRA7A -- for the most part the definitions include the range for the daTRs in addition to flanking TRs
# BOTH:                     283,664                   20: SCA37, NIID, OPML1, DRPLA, SCA2, BSS, FECD, DM1, FAME2, FRA2A, AD, FAME7, FAME3, SCA17, OPDM1, ALS (C9orf72), HSAN8, SBMA, FXS, FRAXE -- the majority of these are defined multiple ways and 1 is the exact definition and the other is the whole region or slightly overlapping
# ISOLATED_only:            3,803,981                 20: HMNR7, FRA11B, FRA12A, SCA3, ALS, FAME6, HDL2, SCA6, OPDM2, PSACH/MED, SPD, GAD, EPM1, SCA10, FAME4, CANVAS, SCA12, FAME1, DMD, XDP
# OVERLAP_no_VC:            312,183                   9: HPE, TOF, BPES, CCHS, SCA1, CCD, HFGS, XLID, XH -- these are all GCN repeats that have subregions defined in addition to the daTR interval (excption: CCD doesn't have an exact definition here!)-- we can add these back individually, I don't want to get into overlapping regions
# BOOKEND_no_VC:            31,218                    3: OPMD, SCA7, FA, 
# TOTAL                     4,542,598                 60 -- all accounted for
# NOTE this total is ~200 regions less than the overall total, may be due to some duplicates in original files
# NOTE this will include some regions that are book-end to be in our overlapping set but they are independent, we would need to revise this to use intersect to report overlapping regions if 
# we want to be more stringent about this -- I am curious how many of the overlapping regions are actually overlapping and how many are just book-ended

# for those still in both and overlapping, can we establish which annotation to keep based on either the TRID or the overlap?
# intersect to get which overlapping regions are subsets of one another
bedtools intersect -C -F 0.90 -a both.bed -b both.bed | sort -k5 -n | awk '{print $NF}' | uniq -c # 199,902 regions are subsets of another region (90% overlap) 
bedtools intersect -C -F 0.90 -a both.bed -b both.bed | awk '$NF > 1 {print}' | bedtools intersect -a ../daTR_coords.sorted.bed -b - # 11/20 of the daTRs in the BOTH set are a larger region that includes another annotation
bedtools intersect -C -F 0.90 -a both.bed -b both.bed | awk '$NF == 1 {print}' | bedtools intersect -a ../daTR_coords.sorted.bed -b - # 18/20 

# If we still want to avoid VCs, it appears like many of them are still genotyped as individual units, so we can take the simplified version of them if there are overlapping annotations
# if we do this how many of the daTRs are preserved? how many overlapping regions remain?
BOTH_vc=${OUT_DIR}/both_vc.bed
BOTH_no_VC=${OUT_DIR}/both_no_vc.bed
BOTH_no_VC_no_overlap=${OUT_DIR}/both_no_vc_no_overlap.bed
BOTH_no_VC_overlap=${OUT_DIR}/both_no_vc_overlap.bed
grep -v --no-group-separator -f /data/projects/nanopore/RepeatExpansion/coordinates/variation_clusters_v1.hg38.TRGT.bed $BOTH > $BOTH_no_VC # this drops 130k entries and leaves 19/20 daTRs in the set
grep --no-group-separator -f /data/projects/nanopore/RepeatExpansion/coordinates/variation_clusters_v1.hg38.TRGT.bed $BOTH > $BOTH_vc
# this is a good way to simplify the annotations and keep the daTRs in the set -- lets keep these separate for now since they do have overlap with VCs so they may be messier to deal with
# but they are the simpiler version of the annotations

# how many of these still overlap each other? we would expect higher than regions that don't overlap with VCs since VCs are just overlapping TRs
bedtools intersect -C -F 0.90 -a $BOTH_no_VC -b $BOTH_no_VC | sort -k5 -n | awk '{print $NF}' | uniq -c  # 151k of these don't overlap another region by 90% or more, OPDM1 has 2 annotations where the correct one has a large overlap with another annotation
bedtools intersect -C -a $BOTH_no_VC -b $BOTH_no_VC | sort -k5 -n | awk '{print $NF}' | uniq -c #125k of these don't overlap another region at all, OPDM1 and SCA2 have 2 overlapping annotations
# it may be alright to include these if we aren't includinf the VCs they overlapped with -- we would still be dropping the 2 daTRs though unless manually added back like the others with overlaps and 
# in VCs
# Of this set, only 120 of them are currently included in our isolated set so adding these would double the daTRs in our set and add 125k regions to the isolated set
bedtools intersect -C -a $BOTH_no_VC -b $BOTH_no_VC | awk '$NF==1 {print $1 "\t" $2 "\t" $3 "\t" $4}' > $BOTH_no_VC_no_overlap # 125k regions that originally overlapped VCs but no other non-VC annotation -- 18 daTRs
bedtools intersect -C -a $BOTH_no_VC -b $BOTH_no_VC | awk '$NF > 1 {print $1 "\t" $2 "\t" $3 "\t" $4}' > $BOTH_no_VC_overlap

# Looking into these sets more closely and how they genotyped I have the following take aways:
# 1. The VCs are going to be very annoying to deal with if we want individual TRs -- they are genotyped as a single unit and we'll have to subset each allele according to the MS
#    There is the additional complication that some of them are defined by multiple motifs but the motif structure isn't stable or some of the motifs never come up
#    This makes this group incredibly heterogeneous and difficult to parse in our current frame work
#    I think it may make the most sense to drop these regions but parse the daTRs (some of which I don't know why they included them in VC when they are only genotyping with 
#    respect to a single motif so they are just isolated TRs being labeled as VCs). There are other daTRs that are less annoying to parse but only defined as VCs (HD, DM2) 
#    but this would be infeasible to do for all of them
# 2. The overlapping regions are mostly overlapping with VCs and are likely to be the same as the VCs -- we can drop these and keep the VCs if we want to keep the VCs and keep the non-overlapping regions
#    regions that remain when we get rid of the VC definitions -- this recovers 18/20 daTRs in this group along with 125k regions (including homopolymers and VNTRS)
# 3. The overlapping regions that don't overlap with VCs are actually overlapping annotations with the exception of ~5k that overlap completely with another annotation
#    this set includes *ALL* 9 of the daTRs in this set which are all GCN repeats that the catalog includes other definitions for -- we can add these back in individually
#    Overall, I don't think we should include overlapping regions in our set because I can't guarantee that I am not double counting their lengths because of poor genotyping in these regions

# overall daTR groups broken down by how we need to address them and if we want to include the sets (*SET* designated the superset, **SET** designates there are multiple levels)
# include key: + include all, * include subset, ** include daTRs, x exclude all

# FILE:                       TOTAL        PERCENT    INCLUDE   daTRs  ACTION
# *VC_only*:                  142,752      3.1%        *     8: SCA8, SCA27B, SCA31, RCPS, SCA36, DM2, HD, FRA7A -- drop these, they are too complicated to parse
# VC_only_multi_motif_STR:    38,463       0.85%       **      
# VC_only_single_motif_STR:   56,605       1.25%       +
# **BOTH**:                   283,664      6.2%        *     20: SCA37, NIID, OPML1, DRPLA, SCA2, BSS, FECD, DM1, FAME2, FRA2A, AD, FAME7, FAME3, SCA17, OPDM1, ALS (C9orf72), HSAN8, SBMA, FXS, FRAXE
# BOTH_VC:                    130,320      2.9%        **    20: AD, ALS.FTD, BSS, DM1, DRPLA, FAME2,FAME3, FAME7, FECD, FRA2A, FRAXE, FXS, HSAN8, NIID, OPDM1, OPML1, SCA17, SCA2, SCA37, SBMA
# *BOTH_no_VC*                153,344      3.4%        *     19: ALS.FTD, BSS, DM1, DRPLA, FAME2,FAME3, FAME7, FECD, FRA2A, FRAXE, FXS, HSAN8, NIID, OPDM1, OPML1, SCA17, SCA2, SCA37, SBMA 
# BOTH_no_VC_no_overlap:      125,955      2.8%        +     17: ALS.FTD, BSS, DM1, DRPLA, FAME2, FAME3, FAME7, FECD, FRA2A, FRAXE, FXS, HSAN8, NIID, OPML1, SBMA, SCA17, SCA37
# BOTH_no_VC_overlap:         27,389       0.6%        **     2: OPDM1, SCA2
# OVERLAP_no_VC:              276,280      6.1%        **     9: HPE, TOF, BPES, CCHS, SCA1, CCD, HFGS, XLID, XH
# BOOKEND_no_VC:              35,903       0.79%       +      3: OPMD, SCA7, FA
# BOOKEND_VC:                 18           0.0%        *      0
# ISOLATED_only:              3,803,981    83.7%       +      20: HMNR7, FRA11B, FRA12A, SCA3, ALS, FAME6, HDL2, SCA6, OPDM2, PSACH/MED, SPD, GAD, EPM1, SCA10, FAME4, CANVAS, SCA12, FAME1, DMD, XDP

# TOTAL INCLUDED:             4,022,455    88.5%             60 -- all accounted for    
# TOTAL EXCLUDED:             520,097      11.4%             0 -- all accounted for
# TOTAL (independent):        4,542,598    100%              60 -- all accounted for

# get the daTR entries so we can get their MEI associations
daTRs=/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/daTR_coords.sorted.bed
bedtools intersect -a $ALL_regions -b $daTRs -wa -wb > $OUT_DIR/daTRs_inputs.bed





