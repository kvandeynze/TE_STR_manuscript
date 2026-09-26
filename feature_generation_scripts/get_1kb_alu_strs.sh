# !/bin/bash
# this script runs bedtools closest to find the closest Alu element and keeps STRs that are within 1kb of an Alu element

# TODO -- should we be limiting this to only STRs that are within 1kb of a full length Alu element? Do both for now, I have the full length Alu subset already generated
# with where their polyA is in the annotation from when I was querying subsequences so I can start there and add in truncated if I want to when I've
# explored how to deal with those -- in theory they should be easier to deal with since if they don't have a polyA the end of their annotation is the end of the Alu element
# however, not all truncated Alus don't have the polyA, they could be 5' truncated and still have the polyA, so I think we should figure out how that affect the 
# relative STR start positions before we even include them here

# first go -- only full length Alus
FULL_LENGTH_ALUS=/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/alu_lineage_data/alu_sets/all_full_length_alus.bed
# This STR set includes all STRs before filtering for sample QC or satellites ect so it should be comparable to the homopolymers and bookend set that we don't apply that to
STRS=/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/all_mei_non_mei_and_daTR_labels_STRs.bed
OUT=/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/random_forest_features/alu_str_features_06192025/full_alu_dists_1kb.bed

# sort both files to ensure bedtools closest works correctly
echo "Sorting $FULL_LENGTH_ALUS and $STRS"

bedtools sort -i $FULL_LENGTH_ALUS > /nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/alu_lineage_data/alu_sets/all_full_length_alus.sorted.bed
bedtools sort -i $STRS > /nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/all_mei_non_mei_and_daTR_labels_STRs.sorted.bed

FULL_LENGTH_ALUS=/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/alu_lineage_data/alu_sets/all_full_length_alus.sorted.bed
STRS=/nfs/boylelab_turbo/kvandeyn/MEI_STR/data/TR_catalog_HPRC/all_mei_non_mei_and_daTR_labels_STRs.sorted.bed

# get the closest Alu element to each STR
echo "Finding closest Alu element to each STR in $STRS"

bedtools closest -a $STRS -b $FULL_LENGTH_ALUS -d | awk -v OFS="\t" '{if ($NF <= 1000) print $0}' > $OUT
echo "Saved closest Alu elements within 1kb of STRs to $OUT"
head -n 10 $OUT