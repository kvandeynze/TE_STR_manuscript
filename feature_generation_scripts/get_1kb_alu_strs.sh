# !/bin/bash
# this script runs bedtools closest to find the closest Alu element and keeps STRs that are within 1kb of an Alu element

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