# IdentifyOffTargetSiteSequences.py
#
# 2015-10-05 Replaced swalign with regex matching.
# 2017-05-31 Replaced nwalign with a explicit search of realignments which uses regex.search.
# 2017-06-03 Output the best offtarget sequences with and/ot without bulges, if any.
# 2025-05-09 Modified alignment to use BioPython PairwiseAligner with precise coordinate mapping.

from __future__ import print_function

import argparse
import collections
import numpy
import os
import string
import operator
import pyfaidx
import re
import subprocess
import logging
from Levenshtein import distance
import pandas as pd
import dill
import sys
from Bio.Align import PairwiseAligner

logger = logging.getLogger("root")


# chromosomePosition defines a class to keep track of the positions.
class chromosomePosition:
    def __init__(self, reference_genome):
        self.chromosome_dict = {}
        self.chromosome_barcode_dict = {}
        self.position_summary = []
        self.index_stack = {}  # we keep track of the values by index here
        self.genome = pyfaidx.Fasta(reference_genome)

    def addPositionBarcode(self, chromosome, position, strand, barcode, primer, count):
        # Create the chromosome keyValue if it doesn't exist
        if chromosome not in self.chromosome_barcode_dict:
            self.chromosome_barcode_dict[chromosome] = {}
        # Increment the position on that chromosome if it exists, otherwise initialize it with 1
        if position not in self.chromosome_barcode_dict[chromosome]:
            self.chromosome_barcode_dict[chromosome][position] = {}
            self.chromosome_barcode_dict[chromosome][position]["+_total"] = 0
            self.chromosome_barcode_dict[chromosome][position]["+primer1_total"] = 0
            self.chromosome_barcode_dict[chromosome][position]["+primer2_total"] = 0
            self.chromosome_barcode_dict[chromosome][position][
                "+primer1_mispriming_total"
            ] = 0
            self.chromosome_barcode_dict[chromosome][position][
                "+primer2_mispriming_total"
            ] = 0
            self.chromosome_barcode_dict[chromosome][position]["+nomatch_total"] = 0

            self.chromosome_barcode_dict[chromosome][position]["-_total"] = 0
            self.chromosome_barcode_dict[chromosome][position]["-primer1_total"] = 0
            self.chromosome_barcode_dict[chromosome][position]["-primer2_total"] = 0
            self.chromosome_barcode_dict[chromosome][position][
                "-primer1_mispriming_total"
            ] = 0
            self.chromosome_barcode_dict[chromosome][position][
                "-primer2_mispriming_total"
            ] = 0
            self.chromosome_barcode_dict[chromosome][position]["-nomatch_total"] = 0

            self.chromosome_barcode_dict[chromosome][position]["+"] = (
                collections.Counter()
            )
            self.chromosome_barcode_dict[chromosome][position]["+primer1"] = (
                collections.Counter()
            )
            self.chromosome_barcode_dict[chromosome][position]["+primer2"] = (
                collections.Counter()
            )
            self.chromosome_barcode_dict[chromosome][position][
                "+primer1_mispriming"
            ] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position][
                "+primer2_mispriming"
            ] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position]["+nomatch"] = (
                collections.Counter()
            )

            self.chromosome_barcode_dict[chromosome][position]["-"] = (
                collections.Counter()
            )
            self.chromosome_barcode_dict[chromosome][position]["-primer1"] = (
                collections.Counter()
            )
            self.chromosome_barcode_dict[chromosome][position]["-primer2"] = (
                collections.Counter()
            )
            self.chromosome_barcode_dict[chromosome][position][
                "-primer1_mispriming"
            ] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position][
                "-primer2_mispriming"
            ] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position]["-nomatch"] = (
                collections.Counter()
            )

        self.chromosome_barcode_dict[chromosome][position][strand][barcode] += count
        self.chromosome_barcode_dict[chromosome][position][strand + primer][
            barcode
        ] += count
        self.chromosome_barcode_dict[chromosome][position][
            strand + primer + "_total"
        ] += count
        self.chromosome_barcode_dict[chromosome][position][strand + "_total"] += count

    def getSequence(self, genome, chromosome, start, end, strand="+"):
        if strand == "+":
            seq = self.genome[chromosome][int(start) : int(end)]
        elif strand == "-":
            seq = self.genome[chromosome][int(start) : int(end)].reverse.complement
        return seq

    # Generates a summary of the barcodes by position
    def SummarizeBarcodePositions(self):
        self.barcode_position_summary = [
            [
                chromosome,
                position,
                len(self.chromosome_barcode_dict[chromosome][position]["+"]),
                len(self.chromosome_barcode_dict[chromosome][position]["-"]),
                self.chromosome_barcode_dict[chromosome][position]["+_total"],
                self.chromosome_barcode_dict[chromosome][position]["-_total"],
                len(self.chromosome_barcode_dict[chromosome][position]["+primer1"]),
                len(self.chromosome_barcode_dict[chromosome][position]["+primer2"]),
                len(self.chromosome_barcode_dict[chromosome][position]["-primer1"]),
                len(self.chromosome_barcode_dict[chromosome][position]["-primer2"]),
                len(
                    self.chromosome_barcode_dict[chromosome][position][
                        "+primer1_mispriming"
                    ]
                ),
                len(
                    self.chromosome_barcode_dict[chromosome][position][
                        "+primer2_mispriming"
                    ]
                ),
                len(
                    self.chromosome_barcode_dict[chromosome][position][
                        "-primer1_mispriming"
                    ]
                ),
                len(
                    self.chromosome_barcode_dict[chromosome][position][
                        "-primer2_mispriming"
                    ]
                ),
            ]
            for chromosome in sorted(self.chromosome_barcode_dict)
            for position in sorted(self.chromosome_barcode_dict[chromosome])
        ]
        self.chr_dataframe_dict = {}
        for c in self.chromosome_barcode_dict:
            self.chr_dataframe_dict[c] = pd.DataFrame(
                self.chromosome_barcode_dict[c].keys()
            )
        return self.barcode_position_summary

    # Summarizes the chromosome, positions within a 10 bp window
    def SummarizeBarcodeIndex(self, windowsize):
        last_chromosome, last_position, window_index = 0, 0, 0
        index_summary = []
        for (
            chromosome,
            position,
            barcode_plus_count,
            barcode_minus_count,
            total_plus_count,
            total_minus_count,
            plus_primer1_count,
            plus_primer2_count,
            minus_primer1_count,
            minus_primer2_count,
            plus_primer1_mispriming_count,
            plus_primer2_mispriming_count,
            minus_primer1_mispriming_count,
            minus_primer2_mispriming_count,
        ) in self.barcode_position_summary:
            if chromosome != last_chromosome or abs(position - last_position) > 10:
                window_index += 1  # new index
            last_chromosome, last_position = chromosome, position
            if window_index not in self.index_stack:
                self.index_stack[window_index] = []
            self.index_stack[window_index].append(
                [
                    chromosome,
                    int(position),
                    int(barcode_plus_count),
                    int(barcode_minus_count),
                    int(barcode_plus_count) + int(barcode_minus_count),
                    int(total_plus_count),
                    int(total_minus_count),
                    int(total_plus_count) + int(total_minus_count),
                    int(plus_primer1_count),
                    int(plus_primer2_count),
                    int(minus_primer1_count),
                    int(minus_primer2_count),
                    int(plus_primer1_mispriming_count),
                    int(plus_primer2_mispriming_count),
                    int(minus_primer1_mispriming_count),
                    int(minus_primer2_mispriming_count),
                ]
            )
        for index in self.index_stack:
            sorted_list = sorted(
                self.index_stack[index], key=operator.itemgetter(4)
            )  # sort by barcode_count_total
            (
                chromosome_list,
                position_list,
                barcode_plus_count_list,
                barcode_minus_count_list,
                barcode_sum_list,
                total_plus_count_list,
                total_minus_count_list,
                total_sum_list,
                plus_primer1_list,
                plus_primer2_list,
                minus_primer1_list,
                minus_primer2_list,
                plus_primer1_mispriming_list,
                plus_primer2_mispriming_list,
                minus_primer1_mispriming_list,
                minus_primer2_mispriming_list,
            ) = zip(*sorted_list)
            barcode_plus = sum(barcode_plus_count_list)
            barcode_minus = sum(barcode_minus_count_list)
            total_plus = sum(total_plus_count_list)
            total_minus = sum(total_minus_count_list)
            plus_primer1 = sum(plus_primer1_list)
            plus_primer2 = sum(plus_primer2_list)
            minus_primer1 = sum(minus_primer1_list)
            minus_primer2 = sum(minus_primer2_list)
            plus_primer1_mispriming = sum(plus_primer1_mispriming_list)
            plus_primer2_mispriming = sum(plus_primer2_mispriming_list)
            minus_primer1_mispriming = sum(minus_primer1_mispriming_list)
            minus_primer2_mispriming = sum(minus_primer2_mispriming_list)
            position_std = numpy.std(position_list)
            min_position = min(position_list)
            max_position = max(position_list)
            barcode_sum = barcode_plus + barcode_minus
            barcode_geometric_mean = (barcode_plus * barcode_minus) ** 0.5
            total_sum = total_plus + total_minus
            total_geometric_mean = (total_plus * total_minus) ** 0.5
            primer1 = plus_primer1 + minus_primer1
            primer2 = plus_primer2 + minus_primer2
            primer1_mispriming = plus_primer1_mispriming + minus_primer1_mispriming
            primer2_mispriming = plus_primer2_mispriming + minus_primer2_mispriming
            primer_geometric_mean = (primer1 * primer2) ** 0.5
            most_frequent_chromosome = sorted_list[-1][0]
            most_frequent_position = sorted_list[-1][1]
            BED_format_chromosome = most_frequent_chromosome
            # if "chr" in most_frequent_chromosome:
            # BED_format_chromosome = most_frequent_chromosome
            # else:
            # BED_format_chromosome = "chr" + most_frequent_chromosome
            # BED_name = BED_format_chromosome + "_" + str(most_frequent_position) + "_" + str(barcode_sum)
            BED_name = f"{BED_format_chromosome}:{min_position}-{max_position}"
            offtarget_sequence = self.getSequence(
                self.genome,
                most_frequent_chromosome,
                most_frequent_position - windowsize,
                most_frequent_position + windowsize,
            )

            summary_list = [
                str(x)
                for x in [
                    index,
                    most_frequent_chromosome,
                    most_frequent_position,
                    offtarget_sequence,  # pick most frequently occurring chromosome and position
                    BED_format_chromosome,
                    min_position,
                    max_position,
                    BED_name,
                    barcode_plus,
                    barcode_minus,
                    barcode_sum,
                    barcode_geometric_mean,
                    total_plus,
                    total_minus,
                    total_sum,
                    total_geometric_mean,
                    primer1,
                    primer2,
                    primer1_mispriming,
                    primer2_mispriming,
                    primer_geometric_mean,
                    position_std,
                ]
            ]

            if barcode_geometric_mean > 0 or primer_geometric_mean > 0:
                index_summary.append(summary_list)
        return index_summary  # WindowIndex, Chromosome, Position, Plus.mi, Minus.mi,
        # BidirectionalArithmeticMean.mi, BidirectionalGeometricMean.mi,
        # Plus, Minus,
        # BidirectionalArithmeticMean, BidirectionalGeometricMean,


"""
Replaces deprecated Bio.pairwise2 with Bio.Align.PairwiseAligner.
Uses glocal alignment (Needleman-Wunsch) with post-processing for Levenshtein distance.
Updated to align Target vs Window and Target vs RC(Window) to maintain strand correctness.
"""


def alignSequences(
    targetsite_sequence, window_sequence, max_score=7, PAM="NGG", max_bulges=2
):
    targetsite_sequence = targetsite_sequence.upper()
    window_sequence = window_sequence.upper()

    # Do NOT reverse complement the target.
    # We look for the Target (5'-3') in the Window (+) and RC-Window (-).

    rc_window_sequence = reverseComplement(window_sequence)

    # Configure PairwiseAligner for Finding Best Anchor (High scores preferred)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2
    # Ensure Global alignment on Query (no skipping allowed)
    # Note: query_end_gap_score defaults to same as gap_score in global mode,
    # which implies penalty. This is correct.

    # Free overhangs for Window (Subject) to allow local matching within window
    aligner.target_end_gap_score = 0.0

    candidates = []

    def process_alignments(query_seq, subject_seq, strand):
        # query_seq = Target Site, subject_seq = Window (or RC Window)
        alignments = aligner.align(query_seq, subject_seq)
        hits = []
        for aln in alignments:
            # aln[0] is aligned Query
            # aln[1] is aligned Subject
            seqA_full = str(aln[0])
            seqB_full = str(aln[1])

            # Identify the core alignment by trimming free end gaps from SeqA (Query)
            match = re.search(r"[^-].*[^-]", seqA_full)
            if not match:
                continue

            start_idx = match.start()
            end_idx = match.end()

            # The 'core' is defined by where the Query aligns.
            coreA = seqA_full[start_idx:end_idx]
            coreB = seqB_full[start_idx:end_idx]

            # Calculate Levenshtein Distance
            subs = 0
            ins = coreB.count("-")  # Gaps in Window (Deletion in Genome/RNA Bulge)
            dels = coreA.count("-")  # Gaps in Target (Insertion in Genome/DNA Bulge)

            for a, b in zip(coreA, coreB):
                if a != "-" and b != "-" and a != b:
                    subs += 1

            total_bulges = ins + dels
            lev_distance = subs + total_bulges

            if lev_distance > max_score:
                continue

            if total_bulges > max_bulges:
                continue

            # Determine start/end in Subject based on the string indices.
            # We map start_idx and end_idx (from the full alignment string)
            # back to the indices of the unaligned Subject string.
            # This is done by counting non-gap characters in seqB up to those points.

            # Count bases in seqB before the start of the core match
            start = len(seqB_full[:start_idx].replace("-", ""))

            # Count bases in seqB before the end of the core match
            end = len(seqB_full[:end_idx].replace("-", ""))

            # Note: end - start will exactly equal the number of bases in coreB (excluding gaps).
            # If coreB contains gaps (RNA bulge), the genomic length will be smaller than Target Length.

            hits.append(
                {
                    "seqA": coreA,  # Aligned Target
                    "seqB": coreB,  # Aligned Window segment (RC if strand -)
                    "score": int(lev_distance),
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "subs": subs,
                    "ins": ins,
                    "del": dels,
                    "bulges": total_bulges,
                }
            )
        return hits

    candidates.extend(process_alignments(targetsite_sequence, window_sequence, "+"))
    candidates.extend(process_alignments(targetsite_sequence, rc_window_sequence, "-"))

    # Initialize return structure
    ret = ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "none"]

    if not candidates:
        return ret

    # Sort candidates: Primary key = Score (Low is better), Secondary = Bulges (Low is better)
    best_hit = sorted(candidates, key=lambda x: (x["score"], x["bulges"]))[0]

    if best_hit["bulges"] == 0:
        ret[0] = best_hit["seqB"]
        ret[1] = best_hit["subs"]
        ret[2] = best_hit["strand"]
        ret[3] = best_hit["start"]
        ret[4] = best_hit["end"]
    else:
        ret[5] = best_hit["seqB"]
        ret[6] = len(best_hit["seqB"])
        ret[7] = best_hit["score"]
        ret[8] = best_hit["subs"]
        ret[9] = best_hit["ins"]
        ret[10] = best_hit["del"]
        ret[11] = best_hit["strand"]
        ret[12] = best_hit["start"]
        ret[13] = best_hit["end"]
        ret[14] = best_hit["seqA"]

    return ret


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def is_control(chr, pos, ref_chr, ref_start, ref_end):
    if chr != ref_chr:
        return False
    if ref_start <= pos <= ref_end:
        return True
    return False


"""
annotation is in the format:
"""


def analyze(
    bam_filename,
    reference_genome,
    outfile,
    annotations,
    windowsize,
    max_score,
    max_bulges,
    control_primer,
    samtools="samtools",
    myDict=None,
):
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    temp = open(outfile + ".primer.tsv", "w")
    tl_filter = open(outfile + ".tl_filter.tsv", "w")
    logger.info("Processing BAM file %s", bam_filename)
    proc = subprocess.Popen(
        [samtools, "view", bam_filename], stdout=subprocess.PIPE, text=True
    )
    file = proc.stdout
    __, filename_tail = os.path.split(bam_filename)
    chromosome_position = chromosomePosition(reference_genome)
    # control_primer_obj = chromosomePosition(reference_genome)
    control_primer_count = 0
    control_counts = 0
    total_dsODN = 0
    control_primer_count_dict = {}  # for debug purposes
    total_dsODN_count_dict = {}  # for debug purposes

    for line in file:
        fields = line.split("\t")
        if len(fields) >= 10:
            # These are strings--need to be cast as ints for comparisons.
            (
                full_read_name,
                sam_flag,
                chromosome,
                position,
                mapq,
                cigar,
                name_of_mate,
                position_of_mate,
                template_length,
                read_sequence,
                read_quality,
            ) = fields[:11]
            if abs(int(template_length)) > 10000:
                print(
                    full_read_name,
                    read_sequence,
                    sam_flag,
                    chromosome,
                    template_length,
                    sep="\t",
                    file=tl_filter,
                )
                continue
            if (
                myDict
                and int(mapq) >= myDict.get("mapq_threshold", 0)
                and int(sam_flag) & 128
                and not int(sam_flag) & 2048
            ):
                # Second read in pair
                barcode, count = parseReadName(full_read_name)
                # print (barcode)
                # print (read_sequence)
                control_read = contain_control_primer(read_sequence, sam_flag, myDict)
                # print (control_read)
                if control_read in control_primer_count_dict:
                    control_primer_count_dict[control_read] += 1
                else:
                    control_primer_count_dict[control_read] = 1
                if control_read != "nomatch":
                    control_primer_count += 1

                # todo, check barcode count
                primer, flag, seq, myDistance, distance2 = assignPrimerstoReads(
                    read_sequence, sam_flag, dsODN_dict=myDict
                )
                if primer != "nomatch":
                    if flag:
                        total_dsODN += 1
                    else:
                        primer = primer + "_mispriming"

                if primer in total_dsODN_count_dict:
                    total_dsODN_count_dict[primer] += 1
                else:
                    total_dsODN_count_dict[primer] = 1
                # print (primer)

                if int(template_length) < 0:  # Reverse read
                    read_position = (
                        int(position_of_mate) + abs(int(template_length)) - 1
                    )
                    strand = "-"
                    chromosome_position.addPositionBarcode(
                        chromosome, read_position, strand, barcode, primer, count
                    )
                elif int(template_length) > 0:  # Forward read
                    read_position = int(position)
                    strand = "+"
                    chromosome_position.addPositionBarcode(
                        chromosome, read_position, strand, barcode, primer, count
                    )
                # if primer == "nomatch":
                print(
                    full_read_name,
                    read_sequence,
                    sam_flag,
                    chromosome,
                    read_position,
                    seq,
                    myDistance,
                    distance2,
                    primer,
                    sep="\t",
                    file=temp,
                )
    proc.stdout.close()
    proc.wait()
    # Generate barcode position summary
    stacked_summary = (
        chromosome_position.SummarizeBarcodePositions()
    )  # this stacked summary is not used
    # print (chromosome_position.chromosome_barcode_dict)
    if control_primer_count == 0:
        control_primer_count = -1
    with open(outfile, "w") as f:
        # Write header
        print(
            "#BED_Chromosome",
            "BED_Min.Position",
            "BED_Max.Position",
            "BED_Name",
            "Filename",
            "WindowIndex",
            "WindowChromosome",
            "Position",
            "WindowSequence",
            "+.mi",
            "-.mi",
            "bi.sum.mi",
            "bi.geometric_mean.mi",
            "+.total",
            "-.total",
            "total.sum",
            "total.geometric_mean",
            "primer1.mi",
            "primer2.mi",
            "primer1_mispriming.mi",
            "primer2_mispriming.mi",
            "primer.geometric_mean",
            "position.stdev",
            "BED_Site_Name",
            "BED_Score",
            "BED_Site_Chromosome",
            "Site_SubstitutionsOnly.Sequence",
            "Site_SubstitutionsOnly.NumSubstitutions",  # 24:25
            "Site_SubstitutionsOnly.Strand",
            "Site_SubstitutionsOnly.Start",
            "Site_SubstitutionsOnly.End",  # 26:28
            "Site_GapsAllowed.Sequence",
            "Site_GapsAllowed.Length",
            "Site_GapsAllowed.Score",  # 29:31
            "Site_GapsAllowed.Substitutions",
            "Site_GapsAllowed.Insertions",
            "Site_GapsAllowed.Deletions",  # 32:34
            "Site_GapsAllowed.Strand",
            "Site_GapsAllowed.Start",
            "Site_GapsAllowed.End",  # 35:37
            "Cell",
            "Targetsite",
            "TargetSequence",
            "RealignedTargetSequence",
            "control_primer_reads",
            "control_counts",
            "normlization_ratio",
            "#pos_500bp",
            "#pos_1kb",
            "#pos_2kb",
            sep="\t",
            file=f,
        )  # 38:41

        # Output summary of each window
        summary = chromosome_position.SummarizeBarcodeIndex(windowsize)
        target_sequence = annotations["Sequence"]
        annotation = [
            annotations["Description"],
            annotations["Targetsite"],
            annotations["Sequence"],
        ]
        output_dict = {}

        for row in summary:
            window_sequence, window_chromosome, window_start, window_end, BED_name = (
                row[3:8]
            )

            non_bulged_target_start_absolute, bulged_target_start_absolute = "", ""
            if target_sequence:
                # Updated to pass max_bulges
                (
                    offtarget_sequence_no_bulge,
                    mismatches,
                    chosen_alignment_strand_m,
                    start_no_bulge,
                    end_no_bulge,
                    bulged_offtarget_sequence,
                    length,
                    distance,
                    substitutions,
                    insertions,
                    deletions,
                    chosen_alignment_strand_b,
                    bulged_start,
                    bulged_end,
                    realigned_target_sequence,
                ) = alignSequences(
                    target_sequence,
                    window_sequence,
                    max_score,
                    myDict["PAM"],
                    max_bulges,
                )

                # print (realigned_target_sequence,target_sequence, window_sequence, max_score)
                BED_score = 1
                BED_chromosome = window_chromosome

                # Logic for absolute coords for mismatch-only alignment
                if chosen_alignment_strand_m == "+":
                    non_bulged_target_start_absolute = (
                        start_no_bulge + int(row[2]) - windowsize
                    )
                    non_bulged_target_end_absolute = (
                        end_no_bulge + int(row[2]) - windowsize
                    )
                elif chosen_alignment_strand_m == "-":
                    non_bulged_target_start_absolute = (
                        int(row[2]) + windowsize - end_no_bulge
                    )
                    non_bulged_target_end_absolute = (
                        int(row[2]) + windowsize - start_no_bulge
                    )
                else:
                    non_bulged_target_start_absolute, non_bulged_target_end_absolute = [
                        ""
                    ] * 2

                # Logic for absolute coords for bulged alignment
                if chosen_alignment_strand_b == "+":
                    bulged_target_start_absolute = (
                        bulged_start + int(row[2]) - windowsize
                    )
                    bulged_target_end_absolute = bulged_end + int(row[2]) - windowsize
                elif chosen_alignment_strand_b == "-":
                    bulged_target_start_absolute = int(row[2]) + windowsize - bulged_end
                    bulged_target_end_absolute = int(row[2]) + windowsize - bulged_start
                else:
                    bulged_target_start_absolute, bulged_target_end_absolute = [""] * 2

                if not (chosen_alignment_strand_m or chosen_alignment_strand_b):
                    BED_chromosome, BED_score, BED_name = [""] * 3
                # print ("BED_name",BED_name)
                # print ("annotation",annotation)
                output_row = (
                    row[4:8]
                    + [filename_tail]
                    + row[0:4]
                    + row[8:]
                    + [
                        str(x)
                        for x in [
                            BED_name,
                            BED_score,
                            BED_chromosome,
                            offtarget_sequence_no_bulge,
                            mismatches,
                            chosen_alignment_strand_m,
                            non_bulged_target_start_absolute,
                            non_bulged_target_end_absolute,
                            bulged_offtarget_sequence,
                            length,
                            distance,
                            substitutions,
                            insertions,
                            deletions,
                            chosen_alignment_strand_b,
                            bulged_target_start_absolute,
                            bulged_target_end_absolute,
                        ]
                    ]
                    + [str(x) for x in annotation]
                    + [realigned_target_sequence]
                )
            else:
                output_row = [
                    str(x)
                    for x in row[4:8]
                    + [filename_tail]
                    + row[0:4]
                    + row[8:]
                    + [""] * 17
                    + annotation
                    + ["none"]
                ]
            # print (output_row)
            if (
                non_bulged_target_start_absolute != ""
                or bulged_target_start_absolute != ""
            ):
                # print ("non_bulged_target_start_absolute",non_bulged_target_start_absolute)
                # print ("bulged_target_start_absolute",bulged_target_start_absolute)
                # print ("non_bulged_target_end_absolute",non_bulged_target_end_absolute)
                # print ("bulged_target_end_absolute",bulged_target_end_absolute)
                output_row_key = "{0}_{1}_{2}".format(
                    window_chromosome,
                    py2min(
                        [non_bulged_target_start_absolute, bulged_target_start_absolute]
                    ),
                    py2max(
                        [non_bulged_target_end_absolute, bulged_target_end_absolute]
                    ),
                )
            else:
                output_row_key = "{0}_{1}_{2}".format(
                    window_chromosome, window_start, window_end
                )

            if output_row_key in output_dict.keys():
                read_count_total = int(output_row[11]) + int(
                    output_dict[output_row_key][11]
                )
                output_dict[output_row_key][11] = str(read_count_total)
            else:
                output_dict[output_row_key] = output_row

        # get control counts:
        for key in sorted(output_dict.keys()):
            current_pos = int(output_dict[key][7])
            chr = output_dict[key][0]
            if (
                myDict
                and "control_coord_chr" in myDict
                and is_control(
                    chr,
                    current_pos,
                    myDict["control_coord_chr"],
                    myDict["control_coord_start"],
                    myDict["control_coord_end"],
                )
            ):
                control_counts = int(float(output_dict[key][11]))
                break
        if control_counts == 0:
            control_counts = -1
        for key in sorted(output_dict.keys()):
            current_pos = int(output_dict[key][7])
            # print (output_dict[key])
            # exit()
            chr = output_dict[key][0]
            # print (chromosome_position.chr_dataframe_dict[chr])
            print(
                *output_dict[key]
                + [
                    str(control_primer_count),
                    str(control_counts),
                    str(float(output_dict[key][11]) / float(control_counts)),
                ]
                + get_num_pos_given_pos(
                    chromosome_position.chr_dataframe_dict[chr], current_pos
                ),
                sep="\t",
                file=f,
            )
    if myDict and myDict.get("save_pickle"):
        save_object(chromosome_position, outfile + ".pkl")


def get_num_pos_given_pos(df, pos):
    # return 500, 1000, and 2000
    out = []
    for i in [500, 1000, 2000]:
        i = int(i / 2)
        out.append(df[0].between(pos - i, pos + i).sum())
    return out


def save_object(obj, out):
    try:
        with open(out, "wb") as f:
            dill.dump(obj, f)
    except Exception as ex:
        print("Error during pickling object (Possibly unsupported):", ex)


def py2min(myList):
    out = [i for i in myList if i != ""]
    return min(out)


def py2max(myList):
    out = [i for i in myList if i != ""]
    return max(out)


def assignPrimerstoReads(read_sequence, sam_flag, dsODN_dict=None):
    if not dsODN_dict:
        return "nomatch", False, "", 100, 100
    # Get 20-nucleotide sequence from beginning or end of sequence depending on orientation
    if int(sam_flag) & 16:
        read_sequence = reverseComplement(read_sequence)
    # i7-
    flag, seq1, myDistance1, distance2 = match_dsODN(
        read_sequence, dsODN_dict["i7-"], dsODN_dict["i7-_match_distance"]
    )
    if flag:
        extend_sequence = read_sequence[
            len(dsODN_dict["i7-"]) : len(dsODN_dict["dsODN_primer"])
        ][: dsODN_dict["i7-_mispriming_length"]]
        correct_sequence = dsODN_dict["dsODN_primer"][
            len(dsODN_dict["i7-"]) : len(dsODN_dict["dsODN_primer"])
        ][: dsODN_dict["i7-_mispriming_length"]]
        # print (extend_sequence,correct_sequence)
        if (
            distance(extend_sequence, correct_sequence)
            >= dsODN_dict["i7-_mispriming_distance"]
        ):
            return "primer1", False, seq1, myDistance1, distance2  # misprimining
        else:
            return "primer1", True, seq1, myDistance1, distance2
    flag, seq2, myDistance2, distance2 = match_dsODN(
        read_sequence, dsODN_dict["i7+"], dsODN_dict["i7+_match_distance"]
    )
    if flag:
        extend_sequence = read_sequence[
            len(dsODN_dict["i7+"]) : len(dsODN_dict["dsODN_primer_revcomp"])
        ][: dsODN_dict["i7+_mispriming_length"]]
        correct_sequence = dsODN_dict["dsODN_primer_revcomp"][
            len(dsODN_dict["i7+"]) : len(dsODN_dict["dsODN_primer_revcomp"])
        ][: dsODN_dict["i7+_mispriming_length"]]
        # print ("i7+",extend_sequence,correct_sequence)
        if (
            distance(extend_sequence, correct_sequence)
            >= dsODN_dict["i7+_mispriming_distance"]
        ):
            return "primer2", False, seq2, myDistance2, distance2  # misprimining
        else:
            return "primer2", True, seq2, myDistance2, distance2
    if myDistance1 < myDistance2:
        return "nomatch", False, seq1, myDistance1, distance2
    return "nomatch", False, seq2, myDistance2, distance2


# i7+_mispriming_length: 12
# i7-_mispriming_length: 5


def match_dsODN(read_sequence, primer, cutoff):
    read_start_sequence = read_sequence[: len(primer)]
    myDistance = distance(read_start_sequence, primer)
    distance2 = myDistance - read_start_sequence.count("N")
    if myDistance <= cutoff:
        return True, read_start_sequence, myDistance, distance2
    else:
        return False, read_start_sequence, myDistance, distance2


def contain_control_primer(read_sequence, sam_flag, myDict=None):
    if not myDict:
        return "nomatch"
    control_primer = myDict["control_primer"]
    if control_primer == "":
        return "nomatch"
    # Get 20-nucleotide sequence from beginning or end of sequence depending on orientation
    if int(sam_flag) & 16:
        read_sequence = reverseComplement(read_sequence)
    myDistance = distance(read_sequence[: len(control_primer)], control_primer)
    if myDistance <= myDict["control_primer_match_distance"]:
        if int(sam_flag) & 16:
            return "control_rev"
        else:
            return "control_fwd"
    return "nomatch"


def loadFileIntoArray(filename):
    with open(filename, "rU") as f:
        keys = f.readline().rstrip("\r\n").split("\t")[1:]
        data = collections.defaultdict(dict)
        for line in f:
            filename, rest = processLine(line)
            line_to_dict = dict(zip(keys, rest))
            data[filename] = line_to_dict
    return data


def parseReadName(read_name):
    # x = read_name.split("_")
    return read_name, 1


def parseReadName2(read_name):
    m = re.search(r"([ACGTN]{8}_[ACGTN]{6}_[ACGTN]{6})_([0-9]*)", read_name)
    if m:
        molecular_index, count = m.group(1), m.group(2)
        return molecular_index, int(count)
    else:
        # print read_name
        return None, None


def processLine(line):
    fields = line.rstrip("\r\n").split("\t")
    filename = fields[0]
    rest = fields[1:]
    return filename, rest


def reverseComplement(sequence):
    if sys.version_info[0] < 3:
        tab = string.maketrans("ACGTacgt", "TGCATGCA")
    else:
        tab = str.maketrans("ACGTacgt", "TGCATGCA")
    return sequence.translate(tab)[::-1]


def main():
    parser = argparse.ArgumentParser(
        description="Identify off-target candidates from Illumina short read sequencing data."
    )
    parser.add_argument("--ref", help="Reference Genome Fasta", required=True)
    parser.add_argument("--samfile", help="SAM file", nargs="*")
    parser.add_argument(
        "--outfile", help="File to output identified sites to.", required=True
    )
    parser.add_argument(
        "--window",
        help="Window around breakpoint to search for off-target",
        type=int,
        default=25,
    )
    parser.add_argument(
        "--control_primer",
        help="control_primer forward sequence",
        type=str,
        default=None,
    )
    parser.add_argument("--max_score", help="Score threshold", type=int, default=7)
    parser.add_argument(
        "--max_bulges", help="Maximum number of bulges allowed", type=int, default=2
    )
    # parser.add_argument('--demo')
    parser.add_argument("--target", default="")

    args = parser.parse_args()

    annotations = {
        "Description": "test description",
        "Targetsite": "dummy targetsite",
        "Sequence": args.target,
    }
    analyze(
        args.samfile[0],
        args.ref,
        args.outfile,
        annotations,
        args.window,
        args.max_score,
        args.max_bulges,
        args.control_primer,
    )


if __name__ == "__main__":
    main()
