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
import logging
from Levenshtein import distance
import pandas as pd
import dill
import sys
from Bio.Align import PairwiseAligner

logger = logging.getLogger("root")


class chromosomePosition:
    def __init__(self, reference_genome):
        self.chromosome_dict = {}
        self.chromosome_barcode_dict = {}
        self.position_summary = []
        self.index_stack = {}
        self.genome = pyfaidx.Fasta(reference_genome)

    def addPositionBarcode(self, chromosome, position, strand, barcode, primer, count):
        if chromosome not in self.chromosome_barcode_dict:
            self.chromosome_barcode_dict[chromosome] = {}
        if position not in self.chromosome_barcode_dict[chromosome]:
            self.chromosome_barcode_dict[chromosome][position] = {}
            for k in [
                "+_total",
                "+primer1_total",
                "+primer2_total",
                "+primer1_mispriming_total",
                "+primer2_mispriming_total",
                "+nomatch_total",
                "-_total",
                "-primer1_total",
                "-primer2_total",
                "-primer1_mispriming_total",
                "-primer2_mispriming_total",
                "-nomatch_total",
            ]:
                self.chromosome_barcode_dict[chromosome][position][k] = 0
            for k in [
                "+",
                "+primer1",
                "+primer2",
                "+primer1_mispriming",
                "+primer2_mispriming",
                "+nomatch",
                "-",
                "-primer1",
                "-primer2",
                "-primer1_mispriming",
                "-primer2_mispriming",
                "-nomatch",
            ]:
                self.chromosome_barcode_dict[chromosome][position][k] = (
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

    def SummarizeBarcodeIndex(self, windowsize):
        last_chromosome, last_position, window_index = 0, 0, 0
        index_summary = []
        for row in self.barcode_position_summary:
            chromosome, position = row[0], row[1]
            if chromosome != last_chromosome or abs(position - last_position) > 10:
                window_index += 1
            last_chromosome, last_position = chromosome, position
            if window_index not in self.index_stack:
                self.index_stack[window_index] = []

            self.index_stack[window_index].append(
                [
                    chromosome,
                    int(position),
                    int(row[2]),
                    int(row[3]),
                    int(row[2]) + int(row[3]),
                    int(row[4]),
                    int(row[5]),
                    int(row[4]) + int(row[5]),
                    int(row[6]),
                    int(row[7]),
                    int(row[8]),
                    int(row[9]),
                    int(row[10]),
                    int(row[11]),
                    int(row[12]),
                    int(row[13]),
                ]
            )

        for index in self.index_stack:
            sorted_list = sorted(self.index_stack[index], key=operator.itemgetter(4))
            data = list(zip(*sorted_list))

            most_frequent_chromosome = sorted_list[-1][0]
            most_frequent_position = sorted_list[-1][1]
            position_list = data[1]
            min_position = min(position_list)
            max_position = max(position_list)
            position_std = numpy.std(position_list)

            barcode_plus = sum(data[2])
            barcode_minus = sum(data[3])
            barcode_sum = barcode_plus + barcode_minus
            barcode_geometric_mean = (barcode_plus * barcode_minus) ** 0.5

            total_plus = sum(data[5])
            total_minus = sum(data[6])
            total_sum = total_plus + total_minus
            total_geometric_mean = (total_plus * total_minus) ** 0.5

            primer1 = sum(data[8]) + sum(data[10])
            primer2 = sum(data[9]) + sum(data[11])
            primer1_mispriming = sum(data[12]) + sum(data[14])
            primer2_mispriming = sum(data[13]) + sum(data[15])
            primer_geometric_mean = (primer1 * primer2) ** 0.5

            BED_format_chromosome = most_frequent_chromosome
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
                    offtarget_sequence,
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
        return index_summary


def alignSequences(
    targetsite_sequence, window_sequence, max_score=7, PAM="NGG", max_bulges=2
):
    targetsite_sequence = targetsite_sequence.upper()
    window_sequence = window_sequence.upper()
    rc_window_sequence = reverseComplement(window_sequence)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2

    # We want full global alignment on the Query (no free end gaps)
    # But we want local alignment on the Window (free end gaps)
    aligner.target_end_gap_score = 0.0

    candidates = []

    def process_alignments(query_seq, subject_seq, strand):
        alignments = aligner.align(query_seq, subject_seq)
        hits = []
        for aln in alignments:
            # We use the full alignment strings.
            # Because query is global, any terminal gaps in seqA are true indels/bulges.
            coreA = str(aln[0])
            coreB = str(aln[1])

            subs = 0
            ins = coreB.count("-")
            dels = coreA.count("-")

            for a, b in zip(coreA, coreB):
                if a != "-" and b != "-" and a != b:
                    subs += 1

            total_bulges = ins + dels
            lev_distance = subs + total_bulges

            if lev_distance > max_score:
                continue
            if total_bulges > max_bulges:
                continue

            # Robust coordinate extraction using aln.aligned
            start, end = 0, 0
            try:
                # aln.aligned[1] contains blocks of indices in the Subject (Window)
                # The first index of the first block is the Start
                # The last index of the last block is the End
                if len(aln.aligned[1]) > 0:
                    start = aln.aligned[1][0][0]
                    end = aln.aligned[1][-1][1]
            except Exception:
                pass

            hits.append(
                {
                    "seqA": coreA,
                    "seqB": coreB,
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

    ret = ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "none"]

    if not candidates:
        return ret

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


def analyze(
    sam_filename,
    reference_genome,
    outfile,
    annotations,
    windowsize,
    max_score,
    max_bulges,
    control_primer,
    myDict=None,
):
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    temp = open(outfile + ".primer.tsv", "w")
    tl_filter = open(outfile + ".tl_filter.tsv", "w")
    logger.info("Processing SAM file %s", sam_filename)
    file = open(sam_filename, "rU")
    __, filename_tail = os.path.split(sam_filename)
    chromosome_position = chromosomePosition(reference_genome)
    control_primer_count = 0
    control_counts = 0
    total_dsODN = 0
    control_primer_count_dict = {}
    total_dsODN_count_dict = {}

    for line in file:
        fields = line.split("\t")
        if len(fields) >= 10:
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
                barcode, count = parseReadName(full_read_name)
                control_read = contain_control_primer(read_sequence, sam_flag, myDict)
                if control_read in control_primer_count_dict:
                    control_primer_count_dict[control_read] += 1
                else:
                    control_primer_count_dict[control_read] = 1
                if control_read != "nomatch":
                    control_primer_count += 1

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

                if int(template_length) < 0:
                    read_position = (
                        int(position_of_mate) + abs(int(template_length)) - 1
                    )
                    strand = "-"
                    chromosome_position.addPositionBarcode(
                        chromosome, read_position, strand, barcode, primer, count
                    )
                elif int(template_length) > 0:
                    read_position = int(position)
                    strand = "+"
                    chromosome_position.addPositionBarcode(
                        chromosome, read_position, strand, barcode, primer, count
                    )
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

    summary = chromosome_position.SummarizeBarcodeIndex(windowsize)
    target_sequence = annotations["Sequence"]
    annotation = [
        annotations["Description"],
        annotations["Targetsite"],
        annotations["Sequence"],
    ]
    output_dict = {}

    for row in summary:
        window_sequence, window_chromosome, window_start, window_end, BED_name = row[
            3:8
        ]
        non_bulged_target_start_absolute, bulged_target_start_absolute = "", ""

        if target_sequence:
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
                target_sequence, window_sequence, max_score, myDict["PAM"], max_bulges
            )

            BED_score = 1
            BED_chromosome = window_chromosome

            if chosen_alignment_strand_m == "+":
                non_bulged_target_start_absolute = (
                    start_no_bulge + int(row[2]) - windowsize
                )
                non_bulged_target_end_absolute = end_no_bulge + int(row[2]) - windowsize
            elif chosen_alignment_strand_m == "-":
                non_bulged_target_start_absolute = (
                    int(row[2]) + windowsize - end_no_bulge
                )
                non_bulged_target_end_absolute = (
                    int(row[2]) + windowsize - start_no_bulge
                )
            else:
                non_bulged_target_start_absolute, non_bulged_target_end_absolute = (
                    "",
                    "",
                )

            if chosen_alignment_strand_b == "+":
                bulged_target_start_absolute = bulged_start + int(row[2]) - windowsize
                bulged_target_end_absolute = bulged_end + int(row[2]) - windowsize
            elif chosen_alignment_strand_b == "-":
                bulged_target_start_absolute = int(row[2]) + windowsize - bulged_end
                bulged_target_end_absolute = int(row[2]) + windowsize - bulged_start
            else:
                bulged_target_start_absolute, bulged_target_end_absolute = "", ""

            if not (chosen_alignment_strand_m or chosen_alignment_strand_b):
                BED_chromosome, BED_score, BED_name = "", "", ""

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

        if non_bulged_target_start_absolute != "" or bulged_target_start_absolute != "":
            output_row_key = "{0}_{1}_{2}".format(
                window_chromosome,
                py2min(
                    [non_bulged_target_start_absolute, bulged_target_start_absolute]
                ),
                py2max([non_bulged_target_end_absolute, bulged_target_end_absolute]),
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

    if control_counts == 0:
        control_counts = -1

    # Calculate control counts first
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

    with open(outfile, "w") as f:
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
            "Site_SubstitutionsOnly.NumSubstitutions",
            "Site_SubstitutionsOnly.Strand",
            "Site_SubstitutionsOnly.Start",
            "Site_SubstitutionsOnly.End",
            "Site_GapsAllowed.Sequence",
            "Site_GapsAllowed.Length",
            "Site_GapsAllowed.Score",
            "Site_GapsAllowed.Substitutions",
            "Site_GapsAllowed.Insertions",
            "Site_GapsAllowed.Deletions",
            "Site_GapsAllowed.Strand",
            "Site_GapsAllowed.Start",
            "Site_GapsAllowed.End",
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
        )
        for key in sorted(output_dict.keys()):
            current_pos = int(output_dict[key][7])
            chr = output_dict[key][0]
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
    if int(sam_flag) & 16:
        read_sequence = reverseComplement(read_sequence)
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
        if (
            distance(extend_sequence, correct_sequence)
            >= dsODN_dict["i7-_mispriming_distance"]
        ):
            return "primer1", False, seq1, myDistance1, distance2
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
        if (
            distance(extend_sequence, correct_sequence)
            >= dsODN_dict["i7+_mispriming_distance"]
        ):
            return "primer2", False, seq2, myDistance2, distance2
        else:
            return "primer2", True, seq2, myDistance2, distance2
    if myDistance1 < myDistance2:
        return "nomatch", False, seq1, myDistance1, distance2
    return "nomatch", False, seq2, myDistance2, distance2


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
    return read_name, 1


def parseReadName2(read_name):
    m = re.search(r"([ACGTN]{8}_[ACGTN]{6}_[ACGTN]{6})_([0-9]*)", read_name)
    if m:
        return m.group(1), int(m.group(2))
    else:
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
