# this is from Matt Danzi -- can probably move this to analysis_utils
import sys
from pathlib import Path

# Make modules in the parent directory importable when this script is run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from analysis_utils import SequenceAttributes
from Levenshtein import distance as edit_distance
import regex as re
from collections import defaultdict


# allele is a string sequence of a single allele call made by TRGT, LongTR, etc
# refMotifs is a list of the reference motifs for this locus according to the catalog
def get_lps_len_motif(allele, refMotifs, maxKmerLen=20):
    """
    Get the length and motif of the longest perfect repeat (LPS) in a given allele and a given set of motifs
    :param allele: str, the allele sequence
    :param refMotifs: list, the reference motifs for this locus according to the catalog
    :return: tuple, (length of the LPS, simplified LPS motif, original motif of the LPS) (modifed so I can have the original motif for purity)
    """
    approvedMotifs={}
    approvedMotifCounts={}
    lpsLength=0
    lpsMotif=''
    spanCovered=0
    done=0
    topLPSLength=0
    topLPSMotif=''
    for k in refMotifs: # this is from the list of reference motifs
        # if any one of the reference motifs has an LPS greater than half the repeat span, that is the LPS
        matches=re.findall('(?:' + k + ')+',allele)
        if len(matches)>0 and len(max(matches))>0.5*len(allele):
            lpsLength = len(max(matches))
            lpsMotif = SequenceAttributes.get_simple_motif(k)
            topLPSMotif = k
            done=1
            break
        else:
            # if not, track the spans covered by this motif
            spanCovered+=sum([len(x) for x in matches])
            if len(matches)>0 and len(max(matches))>topLPSLength:
                topLPSLength=len(max(matches))
                topLPSMotif=k
    if done==0: # this is reached if no reference motif covers more than half the repeat
        # is there enough space left in the repeat for a sequence to be the LPS that isn't already accounted for?
        if (len(allele)-spanCovered)>topLPSLength:
            # do a k-mer search for novel motifs
            if len(allele)<3*maxKmerLen:
                maxKmerLen=int(len(allele)/3)
            for k in range(2,maxKmerLen+1): # changed from 1 to 2, we don't want to return homopolymers
                if done==1:
                    break
                words=re.findall('.'*k,allele)
                words=' '.join(words)
                matches=set(re.findall(r"\b(\w+)\s+\1\b\s\1\b",words)) # only keep k-mers with at least 3 consecutive matches
                for motif in matches:
                    matches=re.findall('(?:' + motif + ')+',allele)
                    spanCovered+=sum([len(x) for x in matches])
                    if len(max(matches))>0.5*len(allele):
                        topLPSLength=len(max(matches))
                        topLPSMotif=motif
                        done=1
                        break
                    elif len(max(matches))>topLPSLength:
                        topLPSLength=len(max(matches))
                        topLPSMotif=motif
        # pick the best LPS
        lpsLength = topLPSLength
        lpsMotif = SequenceAttributes.get_simple_motif(topLPSMotif)

    return lpsLength, lpsMotif, topLPSMotif

def get_lps_len_motif_no_homopolymers(allele, refMotifs, maxKmerLen=20):
    """
    Get the length and motif of the longest perfect repeat (LPS) in a given allele and a given set of motifs
    :param allele: str, the allele sequence
    :param refMotifs: list, the reference motifs for this locus according to the catalog
    :return: tuple, (length of the LPS, simplified LPS motif, original motif of the LPS) (modifed so I can have the original motif for purity)
    """
    approvedMotifs={}
    approvedMotifCounts={}
    lpsLength=0
    lpsMotif=''
    spanCovered=0
    done=0
    topLPSLength=0
    topLPSMotif=''
    for k in refMotifs: # this is from the list of reference motifs
        # if any one of the reference motifs has an LPS greater than half the repeat span, that is the LPS
        matches=re.findall('(?:' + k + ')+',allele)
        if len(matches)>0 and len(max(matches))>0.5*len(allele):
            lpsLength = len(max(matches))
            lpsMotif = SequenceAttributes.get_simple_motif(k)
            topLPSMotif = k
            done=1
            break
        else:
            # if not, track the spans covered by this motif
            spanCovered+=sum([len(x) for x in matches])
            if len(matches)>0 and len(max(matches))>topLPSLength:
                topLPSLength=len(max(matches))
                topLPSMotif=k
    if done==0: # this is reached if no reference motif covers more than half the repeat
        # is there enough space left in the repeat for a sequence to be the LPS that isn't already accounted for?
        if (len(allele)-spanCovered)>topLPSLength:
            # do a k-mer search for novel motifs
            if len(allele)<3*maxKmerLen:
                maxKmerLen=int(len(allele)/3)
            for k in range(2, maxKmerLen + 1):  # changed from 1 to 2, we don't want to return homopolymers
                if done == 1:
                    break
                words = re.findall('.' * k, allele)
                words = ' '.join(words)
                matches = set(re.findall(r"\b(\w+)\s+\1\b\s\1\b", words))  # only keep k-mers with at least 3 consecutive matches
                for motif in matches:
                    # Skip motifs that are composed of only one unique character
                    if len(set(motif)) < 2:
                        continue
                    matches = re.findall('(?:' + motif + ')+', allele)
                    spanCovered += sum([len(x) for x in matches])
                    if len(max(matches)) > 0.5 * len(allele):
                        topLPSLength = len(max(matches))
                        topLPSMotif = motif
                        done = 1
                        break
                    elif len(max(matches)) > topLPSLength:
                        topLPSLength = len(max(matches))
                        topLPSMotif = motif
        # pick the best LPS
        lpsLength = topLPSLength
        lpsMotif = SequenceAttributes.get_simple_motif(topLPSMotif)

    return lpsLength, lpsMotif, topLPSMotif

def fuzzy_segmentation(allele, min_k=2, max_k=12, min_span_ratio=0.1, mismatch_tolerance=4):
    """
    Perform fuzzy segmentation of an allele sequence based on k-mer counting and motif segmentation.
    :param allele: str, the allele sequence
    :param min_k: int, minimum k-mer length
    :param max_k: int, maximum k-mer length
    :param min_span_ratio: float, minimum span ratio for a k-mer to be considered a major k-mer
    :param mismatch_tolerance: int, number of repeat units allowed per mismatch
    :return: list of tuples, each tuple contains (start, end, motif) for each segment
    """

    # Step 1: Perform k-mer counting
    kmer_counts = defaultdict(int)
    for k in range(min_k, max_k + 1):
        for i in range(len(allele) - k + 1):
            kmer = allele[i:i + k]
            kmer_counts[kmer] += 1

    # Step 2: Identify major k-mers
    major_kmers = {}
    for kmer, count in kmer_counts.items():
        span = count * len(kmer)
        if span >= min_span_ratio * len(allele):
            major_kmers[kmer] = span

    # Step 3: Seed fuzzy segmentation
    segments = []
    for motif in major_kmers:
        matches = re.finditer(
            f'(?:{motif}){{1,}}(?:.{{0,{mismatch_tolerance}}}(?:{motif})|(?:.{{0,{mismatch_tolerance}}})){{0,{mismatch_tolerance}}}',
            allele
        )
        for match in matches:
            start, end = match.start(), match.end()
            segments.append((start, end, motif))

    # Step 4: Merge intervals with matching motifs and apply mismatch tolerance
    segments.sort()  # Sort by start position
    merged_segments = []
    for segment in segments:
        if not merged_segments:
            merged_segments.append(segment)
        else:
            prev_start, prev_end, prev_motif = merged_segments[-1]
            curr_start, curr_end, curr_motif = segment

            # Check if motifs match and if the gap between intervals is within mismatch tolerance
            gap = curr_start - prev_end
            if prev_motif == curr_motif and gap <= mismatch_tolerance * len(prev_motif):
                # Merge the intervals
                merged_segments[-1] = (prev_start, curr_end, prev_motif)
            else:
                merged_segments.append(segment)

    # Step 5: Groom boundaries between motifs
    groomed_segments = []
    for i in range(len(merged_segments) - 1):
        start1, end1, motif1 = merged_segments[i]
        start2, end2, motif2 = merged_segments[i + 1]
        if end1 > start2 and motif1.endswith(motif2[:len(motif1)]):
            end1 = start2 + len(motif1) * (start2 // len(motif1))
        groomed_segments.append((start1, end1, motif1))
    groomed_segments.append(merged_segments[-1])  # Add the last segment

    return groomed_segments

def get_purity_lps_motif(lps_motif, allele):
    """
    Get the purity of the LPS motif in a given allele
    This version strictly checkfor the length of pure stretches of the LPS motif
    NOTE: I think this version may do a better job filtering out false positives
    compared to the edit distance version

    :param lps_motif: str, the LPS motif
    :param allele: str, the allele sequence
    :return: float, the purity of the LPS motif in the allele
    """
    matches = re.findall('(?:' + lps_motif + ')+', allele)
    if len(matches) > 0:
        return sum(len(match) for match in matches) / len(allele)
    else:
        return 0.0
def get_purity_lps_edit(allele, lps_motif):
    """
    Get the purity of the LPS motif in a given allele based on edit distance
    :param allele: str, the allele sequence
    :param lps_motif: str, the LPS motif
    :return: float, the purity of the LPS motif in the allele
    """

    if not lps_motif:
        return 0.0

    # Construct the ideal allele made up of repeating lps_motif
    ideal_allele = (lps_motif * (len(allele) // len(lps_motif))) + lps_motif[:len(allele) % len(lps_motif)]

    # Calculate edit distance between the ideal allele and the actual allele
    dist = edit_distance(ideal_allele, allele)

    # Purity is defined as 1 - (edit distance / length of the allele)
    return 1 - (dist / len(allele))