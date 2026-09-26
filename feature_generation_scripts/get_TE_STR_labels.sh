#!/bin/bash

# Script to find closest upstream and downstream repeat elements for each STR
# Reports coordinates, IDs, families, and distances for both upstream and downstream repeats
# Considers LINE, LTR, SINE, and SVA repeat families from RepeatMasker within 25bp of each STR
# Handles overlaps intelligently by classifying relative position

# THIS IS THE SCRIPT I USED TO GENERATE MY FINAL MEI LABELS FOR MY THESIS 

STRS=../QC_scripts/final_updated_labels.bed

REPEATMASKER=/data/projects/nanopore/RepeatExpansion/coordinates/repeatmasker_hg38.fa.out.tab
OUT=.

# Create output directory
mkdir -p $OUT

echo "Extracting LINE, LTR, SINE, and SVA repeat elements from RepeatMasker file..."
# Extract only LINE, LTR, SINE, and SVA repeat elements from RepeatMasker file
# Skip header lines (first 2 lines) and filter for specific families
# Columns: chr(5), start(6), end(7), family(11), ID(15)
# Filter for families starting with LINE, LTR, SINE, or containing SVA
# Ensure all required fields are present and not empty
awk 'NR>2 && NF>=15 && $5!="" && $6!="" && $7!="" && $11!="" && $15!="" && 
     ($11 ~ /^LINE/ || $11 ~ /^LTR/ || $11 ~ /^SINE/ || $11 ~ /SVA/) {
    print $5"\t"$6"\t"$7"\t"$11"\t"$15
}' $REPEATMASKER > ${OUT}/repeat_elements.bed

echo "Number of repeat elements found:"
wc -l ${OUT}/repeat_elements.bed

echo "Repeat family distribution (top 20):"
cut -f4 ${OUT}/repeat_elements.bed | sort | uniq -c | sort -nr | head -20

echo "Sorting files..."
# Validate and clean files before sorting to ensure consistent column numbers
echo "Validating repeat elements file..."
awk 'NF==5' ${OUT}/repeat_elements.bed > ${OUT}/repeat_elements.clean.bed
invalid_repeat_lines=$(awk 'NF!=5' ${OUT}/repeat_elements.bed | wc -l)
echo "Removed $invalid_repeat_lines invalid repeat lines"

echo "Validating STR file..."
awk 'NF==4' $STRS > ${OUT}/strs.clean.bed
invalid_str_lines=$(awk 'NF!=4' $STRS | wc -l)
echo "Removed $invalid_str_lines invalid STR lines"

# Sort both files for bedtools operations
bedtools sort -i ${OUT}/repeat_elements.clean.bed > ${OUT}/repeat_elements.sorted.bed
bedtools sort -i ${OUT}/strs.clean.bed > ${OUT}/strs.sorted.bed

echo "Finding closest upstream repeat elements..."
# Find closest upstream repeat (using -id flag for upstream direction only)
# -D ref reports signed distance relative to reference (negative for upstream)
# First find with overlaps, then without overlaps (-io flag)
bedtools closest -D ref -t all -id -a ${OUT}/strs.sorted.bed -b ${OUT}/repeat_elements.sorted.bed | \
awk '$NF <= 0 && $NF >= -25' > ${OUT}/upstream_repeats_with_overlap.bed

bedtools closest -D ref -t all -id -io -a ${OUT}/strs.sorted.bed -b ${OUT}/repeat_elements.sorted.bed | \
awk '$NF <= 0 && $NF >= -25' > ${OUT}/upstream_repeats_no_overlap.bed

echo "Finding closest downstream repeat elements..."
# Find closest downstream repeat (using -iu flag for downstream direction only)
# -D ref reports signed distance relative to reference (positive for downstream)
# First find with overlaps, then without overlaps (-io flag)
bedtools closest -D ref -t all -iu -a ${OUT}/strs.sorted.bed -b ${OUT}/repeat_elements.sorted.bed | \
awk '$NF >= 0 && $NF <= 25' > ${OUT}/downstream_repeats_with_overlap.bed

bedtools closest -D ref -t all -iu -io -a ${OUT}/strs.sorted.bed -b ${OUT}/repeat_elements.sorted.bed | \
awk '$NF >= 0 && $NF <= 25' > ${OUT}/downstream_repeats_no_overlap.bed

echo "Finding overlapping repeat elements..."
# Find all overlapping repeats to properly classify them
bedtools intersect -wa -wb -a ${OUT}/strs.sorted.bed -b ${OUT}/repeat_elements.sorted.bed > ${OUT}/overlapping_repeats.bed

echo "Number of STRs with upstream repeat within 25bp (with overlaps):"
wc -l ${OUT}/upstream_repeats_with_overlap.bed
echo "Number of STRs with upstream repeat within 25bp (no overlaps):"
wc -l ${OUT}/upstream_repeats_no_overlap.bed

echo "Number of STRs with downstream repeat within 25bp (with overlaps):"
wc -l ${OUT}/downstream_repeats_with_overlap.bed
echo "Number of STRs with downstream repeat within 25bp (no overlaps):"
wc -l ${OUT}/downstream_repeats_no_overlap.bed

echo "Number of STRs with overlapping repeats:"
wc -l ${OUT}/overlapping_repeats.bed

echo "Processing overlapping repeats to determine relative position..."
# For overlapping repeats, determine if they are relatively upstream or downstream
# Based on the center position of the repeat relative to the center of the STR
awk 'BEGIN {OFS="\t"} {
    if (NF >= 9) {
        str_center = ($2 + $3) / 2
        repeat_center = ($6 + $7) / 2
        repeat_start = $6
        repeat_end = $7
        str_start = $2
        str_end = $3
        repeat_family = $8
        repeat_id = $9
        
        # Check if repeat fully encapsulates STR
        if (repeat_start <= str_start && repeat_end >= str_end) {
            # Repeat fully encapsulates STR - report for both upstream and downstream
            print $1, $2, $3, $4, repeat_start, repeat_end, repeat_family, repeat_id, 0, "upstream_overlap_encapsulating"
            print $1, $2, $3, $4, repeat_start, repeat_end, repeat_family, repeat_id, 0, "downstream_overlap_encapsulating"
        } else if (repeat_center < str_center) {
            # Repeat center is upstream of STR center
            print $1, $2, $3, $4, repeat_start, repeat_end, repeat_family, repeat_id, 0, "upstream_overlap"
        } else {
            # Repeat center is downstream of STR center
            print $1, $2, $3, $4, repeat_start, repeat_end, repeat_family, repeat_id, 0, "downstream_overlap"
        }
    }
}' ${OUT}/overlapping_repeats.bed > ${OUT}/overlapping_repeats_classified.bed

echo "Processing upstream repeat data..."
# Process upstream data: combine overlapping upstream + non-overlapping upstream
# For overlapping, use the classified overlaps; for non-overlapping, use the upstream without overlap
awk 'BEGIN {OFS="\t"} $10 ~ /upstream/ {
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, "upstream"
}' ${OUT}/overlapping_repeats_classified.bed > ${OUT}/upstream_repeats_processed.bed

# Add non-overlapping upstream repeats
awk 'BEGIN {OFS="\t"} {
    if (NF >= 10) {
        abs_distance = ($10 < 0) ? -$10 : $10
        print $1, $2, $3, $4, $6, $7, $8, $9, abs_distance, "upstream"
    }
}' ${OUT}/upstream_repeats_no_overlap.bed >> ${OUT}/upstream_repeats_processed.bed

echo "Processing downstream repeat data..."
# Process downstream data: combine overlapping downstream + non-overlapping downstream
awk 'BEGIN {OFS="\t"} $10 ~ /downstream/ {
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, "downstream"
}' ${OUT}/overlapping_repeats_classified.bed > ${OUT}/downstream_repeats_processed.bed

# Add non-overlapping downstream repeats
awk 'BEGIN {OFS="\t"} {
    if (NF >= 10) {
        print $1, $2, $3, $4, $6, $7, $8, $9, $10, "downstream"
    }
}' ${OUT}/downstream_repeats_no_overlap.bed >> ${OUT}/downstream_repeats_processed.bed

echo "Combining upstream and downstream results..."
# Combine both upstream and downstream results
cat ${OUT}/upstream_repeats_processed.bed ${OUT}/downstream_repeats_processed.bed > ${OUT}/flanking_repeats_combined.bed

# Sort by STR coordinates for easier analysis
bedtools sort -i ${OUT}/flanking_repeats_combined.bed > ${OUT}/flanking_repeats_combined.sorted.bed

echo "Creating summary by STR..."
# Create a summary showing both upstream and downstream repeats for each STR
# For each STR, find the closest upstream and downstream repeat
# Handle encapsulating cases specially
# IMPORTANT: Process ALL STRs, not just those with flanking repeats

# First, collect flanking information for STRs that have it
awk 'BEGIN {OFS="\t"}
{
    str_key = $1"_"$2"_"$3"_"$4
    direction = $10
    distance = $9
    repeat_coords = $5":"$6
    repeat_family = $7
    repeat_id = $8
    
    if (direction == "upstream" || direction == "upstream_overlap" || direction == "upstream_overlap_encapsulating") {
        # Keep track of closest upstream (smallest distance)
        if (!(str_key in upstream_dist) || distance < upstream_dist[str_key]) {
            upstream_coords[str_key] = repeat_coords
            upstream_family[str_key] = repeat_family
            upstream_id[str_key] = repeat_id
            upstream_dist[str_key] = distance
        }
        # Special handling for encapsulating - also store as downstream
        if (direction == "upstream_overlap_encapsulating") {
            if (!(str_key in downstream_dist) || distance < downstream_dist[str_key]) {
                downstream_coords[str_key] = repeat_coords
                downstream_family[str_key] = repeat_family
                downstream_id[str_key] = repeat_id
                downstream_dist[str_key] = distance
            }
        }
    } else if (direction == "downstream" || direction == "downstream_overlap" || direction == "downstream_overlap_encapsulating") {
        # Keep track of closest downstream (smallest distance)
        if (!(str_key in downstream_dist) || distance < downstream_dist[str_key]) {
            downstream_coords[str_key] = repeat_coords
            downstream_family[str_key] = repeat_family
            downstream_id[str_key] = repeat_id
            downstream_dist[str_key] = distance
        }
        # Special handling for encapsulating - also store as upstream
        if (direction == "downstream_overlap_encapsulating") {
            if (!(str_key in upstream_dist) || distance < upstream_dist[str_key]) {
                upstream_coords[str_key] = repeat_coords
                upstream_family[str_key] = repeat_family
                upstream_id[str_key] = repeat_id
                upstream_dist[str_key] = distance
            }
        }
    }
}
END {
    # Output flanking data for processing by next step
    for (str_key in upstream_coords) {
        print str_key, upstream_coords[str_key], upstream_family[str_key], upstream_id[str_key], upstream_dist[str_key], "upstream"
    }
    for (str_key in downstream_coords) {
        print str_key, downstream_coords[str_key], downstream_family[str_key], downstream_id[str_key], downstream_dist[str_key], "downstream"
    }
}' ${OUT}/flanking_repeats_combined.sorted.bed > ${OUT}/flanking_data.temp

# Now process ALL STRs from the original file and add flanking information where available
awk 'BEGIN {OFS="\t"}
# Read flanking data first
FNR==NR {
    str_key = $1
    if ($6 == "upstream") {
        upstream_coords[str_key] = $2
        upstream_family[str_key] = $3
        upstream_id[str_key] = $4
        upstream_dist[str_key] = $5
    } else if ($6 == "downstream") {
        downstream_coords[str_key] = $2
        downstream_family[str_key] = $3
        downstream_id[str_key] = $4
        downstream_dist[str_key] = $5
    }
    next
}
# Process all STRs
{
    if (NF >= 4) {
        str_key = $1"_"$2"_"$3"_"$4
        
        up_coords = (upstream_coords[str_key] != "") ? upstream_coords[str_key] : "NA"
        up_family = (upstream_family[str_key] != "") ? upstream_family[str_key] : "NA"
        up_id = (upstream_id[str_key] != "") ? upstream_id[str_key] : "NA"
        up_dist = (upstream_dist[str_key] != "") ? upstream_dist[str_key] : "NA"
        down_coords = (downstream_coords[str_key] != "") ? downstream_coords[str_key] : "NA"
        down_family = (downstream_family[str_key] != "") ? downstream_family[str_key] : "NA"
        down_id = (downstream_id[str_key] != "") ? downstream_id[str_key] : "NA"
        down_dist = (downstream_dist[str_key] != "") ? downstream_dist[str_key] : "NA"
        
        # Classify STRs based on MEI relationship
        classification = "none"
        
        # Case 1: STR has both upstream and downstream MEI information
        if (up_id != "NA" && down_id != "NA") {
            # Internal: Same MEI ID upstream and downstream (fragments of same MEI)
            if (up_id == down_id) {
                # Check if fully overlapped (same coordinates) or between fragments (different coordinates)
                if (up_coords == down_coords && up_dist == "0" && down_dist == "0") {
                    classification = "internal_fully_overlapped"
                } else if (up_coords != down_coords) {
                    classification = "internal_between_fragments"
                } else {
                    classification = "internal_overlapped"
                }
            }
            # External: Different MEI IDs upstream and downstream
            else {
                classification = "external_different_meis"
            }
        }
        # Case 2: STR has only upstream MEI
        else if (up_id != "NA" && down_id == "NA") {
            if (up_dist == "0") {
                classification = "external_upstream_overlap"
            } else {
                classification = "external_upstream_proximity"
            }
        }
        # Case 3: STR has only downstream MEI
        else if (up_id == "NA" && down_id != "NA") {
            if (down_dist == "0") {
                classification = "external_downstream_overlap"
            } else {
                classification = "external_downstream_proximity"
            }
        }
        # Case 4: No MEI within 25bp
        else {
            classification = "none"
        }
        
        print $1, $2, $3, $4, up_coords, up_family, up_id, up_dist, down_coords, down_family, down_id, down_dist, classification
    }
}' ${OUT}/flanking_data.temp ${OUT}/strs.sorted.bed > ${OUT}/flanking_repeats_summary_temp.tsv

# Add header and create final summary file
echo -e "STR_chr\tSTR_start\tSTR_end\tSTR_ID\tupstream_coords\tupstream_family\tupstream_ID\tupstream_distance\tdownstream_coords\tdownstream_family\tdownstream_ID\tdownstream_distance\tSTR_classification" > ${OUT}/flanking_repeats_summary.tsv
cat ${OUT}/flanking_repeats_summary_temp.tsv >> ${OUT}/flanking_repeats_summary.tsv

# Clean up temporary files
rm ${OUT}/flanking_data.temp ${OUT}/flanking_repeats_summary_temp.tsv

echo "Results summary:"
echo "Total STRs processed: $(wc -l < ${OUT}/strs.sorted.bed)"
echo "STRs with overlapping repeats: $(wc -l < ${OUT}/overlapping_repeats.bed)"
echo "STRs with upstream repeat within 25bp: $(tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$5 != "NA"' | wc -l)"
echo "STRs with downstream repeat within 25bp: $(tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$9 != "NA"' | wc -l)"
echo "STRs with either upstream or downstream repeat: $(tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$5 != "NA" || $9 != "NA"' | wc -l)"
echo "STRs with both upstream and downstream repeats: $(tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$5 != "NA" && $9 != "NA"' | wc -l)"
echo "STRs with encapsulating repeats (same repeat for both upstream and downstream): $(tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$7 == $11 && $7 != "NA"' | wc -l)"

echo ""
echo "STR Classification Summary:"
echo "==========================="
echo "INTERNAL STRs (fully enclosed by MEI or between MEI fragments):"
tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$13 ~ /^internal/' | cut -f13 | sort | uniq -c | sort -nr

echo ""
echo "EXTERNAL STRs (partial overlap or proximity to MEI):"
tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$13 ~ /^external/' | cut -f13 | sort | uniq -c | sort -nr

echo ""
echo "Summary by main categories:"
echo "Internal STRs: $(tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$13 ~ /^internal/' | wc -l)"
echo "External STRs: $(tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$13 ~ /^external/' | wc -l)"
echo "No MEI within 25bp: $(tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$13 == "none"' | wc -l)"

echo "Repeat family distribution in flanking results:"
echo "Upstream families (top 15):"
tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$6 != "NA" {print $6}' | sort | uniq -c | sort -nr | head -15

echo "Downstream families (top 15):"
tail -n +2 ${OUT}/flanking_repeats_summary.tsv | awk '$10 != "NA" {print $10}' | sort | uniq -c | sort -nr | head -15

echo ""
echo "Output files created:"
echo "1. ${OUT}/flanking_repeats_combined.sorted.bed - All upstream and downstream repeats with distances"
echo "2. ${OUT}/flanking_repeats_summary.tsv - Summary table with upstream and downstream info per STR"
echo "3. ${OUT}/overlapping_repeats_classified.bed - Overlapping repeats classified by relative position"
echo ""
echo "Column structure for flanking_repeats_summary.tsv:"
echo "1. STR_chr - STR chromosome"
echo "2. STR_start - STR start position"
echo "3. STR_end - STR end position"
echo "4. STR_ID - STR LocusID"
echo "5. upstream_coords - Upstream repeat coordinates (start:end)"
echo "6. upstream_family - Upstream repeat family/class"
echo "7. upstream_ID - Upstream repeat RepeatMasker ID"
echo "8. upstream_distance - Distance to upstream repeat (0 for overlapping)"
echo "9. downstream_coords - Downstream repeat coordinates (start:end)"
echo "10. downstream_family - Downstream repeat family/class"
echo "11. downstream_ID - Downstream repeat RepeatMasker ID"
echo "12. downstream_distance - Distance to downstream repeat (0 for overlapping)"
echo "13. STR_classification - Classification of STR relative to MEI elements"
echo "9. downstream_coords - Downstream repeat coordinates (start:end)"
echo "10. downstream_family - Downstream repeat family/class"
echo "11. downstream_ID - Downstream repeat RepeatMasker ID"
echo "12. downstream_distance - Distance to downstream repeat (0 for overlapping)"

echo ""
echo "=============================================================================="
echo "DETAILED OUTPUT FILE DESCRIPTIONS"
echo "=============================================================================="
echo ""

echo "1. flanking_repeats_summary.tsv (Main Results File)"
echo "   - Tab-separated values with header"
echo "   - One row per STR with flanking repeat information"
echo "   - Columns:"
echo "     1. STR_chr: Chromosome of the STR"
echo "     2. STR_start: Start coordinate of the STR"
echo "     3. STR_end: End coordinate of the STR"
echo "     4. STR_ID: Unique identifier for the STR (LocusID)"
echo "     5. upstream_coords: Coordinates of closest upstream repeat (format: start:end) or 'NA'"
echo "     6. upstream_family: RepeatMasker family/class of upstream repeat or 'NA'"
echo "     7. upstream_ID: RepeatMasker ID of upstream repeat or 'NA'"
echo "     8. upstream_distance: Distance in bp to upstream repeat (0 for overlapping) or 'NA'"
echo "     9. downstream_coords: Coordinates of closest downstream repeat (format: start:end) or 'NA'"
echo "     10. downstream_family: RepeatMasker family/class of downstream repeat or 'NA'"
echo "     11. downstream_ID: RepeatMasker ID of downstream repeat or 'NA'"
echo "     13. STR_classification: Classification of STR relative to MEI elements - one of:"
echo "        INTERNAL classifications (STR fully enclosed by MEI):"
echo "          - 'internal_fully_overlapped': STR completely overlapped by single MEI fragment"
echo "          - 'internal_between_fragments': STR between upstream and downstream fragments of same MEI"
echo "          - 'internal_overlapped': STR overlapped by same MEI upstream and downstream"
echo "        EXTERNAL classifications (STR partially overlapped or near MEI):"
echo "          - 'external_different_meis': STR flanked by different MEI elements"
echo "          - 'external_upstream_overlap': STR overlaps MEI only upstream"
echo "          - 'external_downstream_overlap': STR overlaps MEI only downstream"
echo "          - 'external_upstream_proximity': STR near MEI upstream (within 25bp, no overlap)"
echo "          - 'external_downstream_proximity': STR near MEI downstream (within 25bp, no overlap)"
echo "        OTHER:"
echo "          - 'none': No MEI within 25bp of STR"
echo ""

echo "2. flanking_repeats_combined.sorted.bed (Detailed Results File)"
echo "   - BED-like format, tab-separated, no header"
echo "   - One row per STR-repeat flanking relationship"
echo "   - Sorted by genomic coordinates"
echo "   - Columns:"
echo "     1. STR_chr: Chromosome of the STR"
echo "     2. STR_start: Start coordinate of the STR"
echo "     3. STR_end: End coordinate of the STR"
echo "     4. STR_ID: Unique identifier for the STR (LocusID)"
echo "     5. Repeat_start: Start coordinate of the repeat element"
echo "     6. Repeat_end: End coordinate of the repeat element"
echo "     7. Repeat_family: RepeatMasker family/class of the repeat element"
echo "     8. Repeat_ID: RepeatMasker ID of the repeat element"
echo "     9. distance: Distance in bp between STR and repeat (0 for overlapping)"
echo "     10. direction: Direction relative to STR ('upstream' or 'downstream')"
echo ""

echo "3. overlapping_repeats_classified.bed (Overlap Classification File)"
echo "   - BED-like format, tab-separated, no header"
echo "   - One row per STR-repeat overlap relationship"
echo "   - Shows how overlapping repeats were classified"
echo "   - Columns:"
echo "     1. STR_chr: Chromosome of the STR"
echo "     2. STR_start: Start coordinate of the STR"
echo "     3. STR_end: End coordinate of the STR"
echo "     4. STR_ID: Unique identifier for the STR (LocusID)"
echo "     5. Repeat_start: Start coordinate of the overlapping repeat element"
echo "     6. Repeat_end: End coordinate of the overlapping repeat element"
echo "     7. Repeat_family: RepeatMasker family/class of the overlapping repeat"
echo "     8. Repeat_ID: RepeatMasker ID of the overlapping repeat element"
echo "     9. distance: Always 0 for overlapping elements"
echo "     10. classification: One of:"
echo "        - 'upstream_overlap': Repeat center is upstream of STR center"
echo "        - 'downstream_overlap': Repeat center is downstream of STR center"
echo "        - 'upstream_overlap_encapsulating': Repeat fully contains STR (reported as upstream)"
echo "        - 'downstream_overlap_encapsulating': Repeat fully contains STR (reported as downstream)"
echo ""

echo "4. Intermediate Files (for debugging/analysis):"
echo "   - repeat_elements.bed: All repeat elements extracted from RepeatMasker"
echo "   - repeat_elements.clean.bed: Cleaned repeat elements (consistent 5 columns)"
echo "   - repeat_elements.sorted.bed: Sorted repeat elements for bedtools"
echo "   - strs.clean.bed: Cleaned STR file (consistent 4 columns)"
echo "   - strs.sorted.bed: Sorted STR file for bedtools"
echo "   - upstream_repeats_with_overlap.bed: Upstream repeats including overlaps"
echo "   - upstream_repeats_no_overlap.bed: Upstream repeats excluding overlaps"
echo "   - downstream_repeats_with_overlap.bed: Downstream repeats including overlaps"
echo "   - downstream_repeats_no_overlap.bed: Downstream repeats excluding overlaps"
echo "   - overlapping_repeats.bed: Raw overlapping repeat-STR pairs from bedtools intersect"
echo ""

echo "=============================================================================="
echo "ANALYSIS NOTES:"
echo "=============================================================================="
echo "- Distance cutoff: 25 bp maximum for both upstream and downstream"
echo "- Repeat families: LINE, LTR, SINE, and SVA families only (MEI-related repeats)"
echo "- Overlap handling: Overlapping repeats are classified by relative position"
echo "- Encapsulation: When repeat fully contains STR, same repeat reported for both directions"
echo "- Priority: For multiple repeats at same distance, closest by coordinate is chosen"
echo "- Missing data: 'NA' indicates no repeat found within 25bp in that direction"
echo "- Distance calculation: 0 for overlapping, actual bp distance for non-overlapping"
echo "- Family information: Includes RepeatMasker family/class for each flanking repeat"
echo "=============================================================================="
