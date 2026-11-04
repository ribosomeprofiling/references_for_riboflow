
import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import gzip
from collections import defaultdict
import statistics

# Not the best way of imporing modules.
# we hope to organize our scripts and librareis better later.
script_path = os.path.dirname( os.path.realpath(__file__) )
#upper_path =  os.path.dirname(os.path.dirname(script_path))
lib_path = os.path.join(script_path,"..","..","..", "..")
if not os.path.isdir(lib_path):
    print("Couldn't find the dir path {}. Exiting...".format(lib_path))
    exit(1)

sys.path.insert(0, lib_path)

import ref_lib
from ref_lib.Fasta import FastaFile, FastaEntry
from ref_lib.GTF import GTFfile, GTFEntry, get_gtf_contents_genome

###############################################################################
# Region Names:
UTR5_name          = "UTR5"
UTR5_junction_name = "UTR5_junction"
CDS_name           = "CDS"
UTR3_junction_name = "UTR3_junction"
UTR3_name          = "UTR3"

start_site_name = "start_site"
stop_site_name = "stop_site"

###############################################################################

# TODO:
# Add functions to display statistics

###########################################################

def get_arguments():
    parser = argparse.ArgumentParser(description="""


    """, formatter_class = RawTextHelpFormatter)

    parser.add_argument("-f" ,
                        help = """Gencode fasta File Containing Transcript Sequences""" ,
                        required = True ,
                        type = str)
    parser.add_argument("-g" ,
                        help = """Gencode gtf file. This file contains genome annotations""" ,
                        required = True ,
                        type = str)
    parser.add_argument("-a" ,
                        help = """Appris isoform list""" ,
                        required = True ,
                        type = str)
    parser.add_argument("-o" ,
                        help = """Output base.""" ,
                        required = True ,
                        type = str)

    return parser.parse_args()

###############################################################

def _find_rel_cds_end_coords( exon_list, cds_contents, coord_type, strand ):
    """
    Returns the start or end of the CDS 0-based bed format
    So the start position is included and end position is excluded
    """
    # c is the corresponding_exon_index to the cds chunk

    if coord_type == "start":
        cds_chunk = cds_contents[0]
        c = cds_chunk[2]
        if strand == "+":
            this_sum = cds_chunk[0] - exon_list[c][0]
        else:
            this_sum = exon_list[c][1] - cds_chunk[1]
    else:
        cds_chunk = cds_contents[-1]
        c = cds_chunk[2]
        this_sum = (cds_chunk[1] - cds_chunk[0]) + 1

        if len(cds_contents) == 1:
            if strand == "+":
                this_sum += cds_chunk[0]  - exon_list[c][0]
            else:
                this_sum += exon_list[c][1] - cds_chunk[1]
        """
        Delete this later!!!
        if strand == "+":
            this_sum = (cds_chunk[1] - cds_chunk[0]) + 1
        else:
            this_sum = cds_chunk[0] - exon_list[c][0]
        """
    for e in exon_list[:c]:
        this_sum += (e[1] - e[0]) + 1

    return this_sum

###########################################################################################

def find_rel_utr_cds_positions( genes_dict ):
    """
    For each transcri,pt, the first nucleotide (of the transcript) has position 0
    the coordinates are in BED stye
    so For A,B
    A: start, is INCLUDED
    B: end, is EXCLUDED
    The coordinates have the _bed prefix and they are stored in
    the corresponding transcript dictionary
    """
    #
    # 0-based and end coordinate is excluded
    #
    for g_name, transcripts in genes_dict.items():

        for t_name, t_contents in transcripts.items():
            CDS_contents = t_contents.get("CDS", list())
            if len(CDS_contents) == 0:
                continue

            cds_rel_start = \
                _find_rel_cds_end_coords( t_contents["exons"], CDS_contents, "start", t_contents["strand"]   )
            cds_rel_end = \
                _find_rel_cds_end_coords( t_contents["exons"], CDS_contents, "end", t_contents["strand"] )

            t_length = t_contents["length"]

            t_contents["bed_CDS"] = ( cds_rel_start, cds_rel_end )
            if cds_rel_start > 0:
                t_contents["bed_UTR_5"] = (0, cds_rel_start  )

            if cds_rel_end < t_length:
                t_contents["bed_UTR_3"] = (cds_rel_end, t_length )


#############################################################################################

def read_principal_appris_trascript_list(appris_file, gtf_contents):
    """
    Reads the appris isoform file into a dictionary

    We assume that the files is of the form:
    Gene NAme \t GENE_ID \t TRANSCRIOPT_ID \t APPRIS_CATEGORY \t MANE
    # C1orf112  ENSG00000000460 ENST00000286031 CCDS1285    PRINCIPAL:1  MANE_Select

    The output dictionary is of the form:
    appris_genes[gene_id][transcript_id] = {  "gene_name": gene_name, "category": category }

    Also returns a set of (gene_id, transcript_id) tuples for transcripts
    that are both PRINCIPAL and MANE_Select.
    """

    appris_genes = defaultdict(dict)
    mane_select_transcripts = set()
    myopen=open
    if appris_file.endswith(".gz"):
        myopen = gzip.open

    with myopen(appris_file, "rt") as input_stream:
        for entry in input_stream:
            contents = entry.strip().split("\t")
            if len(contents) < 5:
                continue
            gene_name, gene_id, transcript_id, dummy, category = contents[:5]
            mane = contents[5].strip() if len(contents) > 5 else ""
            if not category.startswith("PRINCIPAL"):
                continue
            level_str = category.split(":")[1]
            principal_level = int(level_str) if level_str.isdigit() else 99
            appris_genes[gene_id][transcript_id] = gtf_contents[gene_id][transcript_id]
            appris_genes[gene_id][transcript_id]["principal_level"] = principal_level
            appris_genes[gene_id][transcript_id]["gene_name"] = gene_name
            if mane == "MANE_Select":
                mane_select_transcripts.add((gene_id, transcript_id))

    return appris_genes, mane_select_transcripts


#############################################################################################

def pick_best_appris_transcripts(appris_genes, mane_select_transcripts):
    """
    For each gene, pick one transcript using the following priority:
    1. PRINCIPAL:1 + MANE_Select
    2. PRINCIPAL:2 + MANE_Select
    3. PRINCIPAL:1
    4. PRINCIPAL:2
    5. MANE_Select (any level)
    6. Longest PRINCIPAL (any level)
    Ties within a tier are broken by longest transcript length.
    """
    best_picks = defaultdict(dict)
    counts = {
        "P1+MANE": 0,
        "P2+MANE": 0,
        "P1": 0,
        "P2": 0,
        "MANE": 0,
        "longest": 0,
    }

    def longest(candidates):
        return max(candidates, key=lambda x: x[1]["length"])

    for gene, transcripts in appris_genes.items():
        t_items = list(transcripts.items())

        is_mane  = lambda t: (gene, t) in mane_select_transcripts
        by_level = lambda lvl: [(t, c) for t, c in t_items if c["principal_level"] == lvl]

        p1_mane = [(t, c) for t, c in t_items if c["principal_level"] == 1 and is_mane(t)]
        p2_mane = [(t, c) for t, c in t_items if c["principal_level"] == 2 and is_mane(t)]
        p1      = by_level(1)
        p2      = by_level(2)
        mane    = [(t, c) for t, c in t_items if is_mane(t)]

        if p1_mane:
            best = longest(p1_mane)
            counts["P1+MANE"] += 1
        elif p2_mane:
            best = longest(p2_mane)
            counts["P2+MANE"] += 1
        elif p1:
            best = longest(p1)
            counts["P1"] += 1
        elif p2:
            best = longest(p2)
            counts["P2"] += 1
        elif mane:
            best = longest(mane)
            counts["MANE"] += 1
        else:
            best = longest(t_items)
            counts["longest"] += 1

        best_picks[gene][best[0]] = best[1]

    print("  Picked by PRINCIPAL:1 + MANE_Select: {}".format(counts["P1+MANE"]))
    print("  Picked by PRINCIPAL:2 + MANE_Select: {}".format(counts["P2+MANE"]))
    print("  Picked by PRINCIPAL:1 only:           {}".format(counts["P1"]))
    print("  Picked by PRINCIPAL:2 only:           {}".format(counts["P2"]))
    print("  Picked by MANE_Select:                {}".format(counts["MANE"]))
    print("  Picked by longest PRINCIPAL:          {}".format(counts["longest"]))
    print("  Total picks: {}".format(len(best_picks)))
    return best_picks


#############################################################################################

def read_readthrough_gene_ids(gtf_file):
    """
    Returns the set of gene ids (version suffix stripped) whose gene entry
    carries the GENCODE tag "readthrough_gene".
    """
    readthrough_genes = set()
    myopen = gzip.open if gtf_file.endswith(".gz") else open

    with myopen(gtf_file, "rt") as input_stream:
        for line in input_stream:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] != "gene":
                continue
            if 'tag "readthrough_gene"' not in fields[8]:
                continue
            gene_id = fields[8].split('gene_id "')[1].split('"')[0]
            readthrough_genes.add(gene_id.split(".")[0])

    return readthrough_genes


def remove_readthrough_genes(transcript_dict, readthrough_gene_ids):
    """
    Drops genes tagged as read-through from the transcript dictionary.
    Read-through transcripts duplicate exon sequence of their parent genes,
    creating multi-mapping ambiguity in the transcriptome reference.
    """
    kept = defaultdict(dict)
    removed_count = 0

    for g, transcripts in transcript_dict.items():
        if g in readthrough_gene_ids:
            removed_count += 1
        else:
            kept[g] = transcripts

    print("  Removed {} read-through genes.".format(removed_count))
    return kept

#############################################################################################

def _merge_intervals(intervals):
    """Merges overlapping/adjacent (start, end) intervals; coordinates inclusive."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def _exons_contained_in(inner_exons, outer_exons):
    """True if every inner interval lies within some outer interval (both merged, sorted)."""
    o = 0
    for start, end in inner_exons:
        while o < len(outer_exons) and outer_exons[o][1] < end:
            o += 1
        if o == len(outer_exons) or outer_exons[o][0] > start:
            return False
    return True


def remove_contained_transcripts(transcript_dict):
    """
    Removes transcripts whose exonic bases all lie within the merged exons of
    another, strictly longer transcript on the same chromosome and strand.
    These duplicate their container's sequence in the reference, so reads from
    them are ambiguous. Transcripts with identical exon structure (equal
    length) are both kept.
    """
    groups = defaultdict(list)
    for g, transcripts in transcript_dict.items():
        for t, c in transcripts.items():
            merged = _merge_intervals(c["exons"])
            groups[(c["chrom"], c["strand"])].append(
                (merged[0][0], merged[-1][1], g, t, merged, c["length"]))

    removed = set()
    for items in groups.values():
        items.sort()
        active = []
        for cur in items:
            active = [a for a in active if a[1] >= cur[0]]
            for other in active:
                for inner, outer in ((cur, other), (other, cur)):
                    if inner[5] < outer[5] and _exons_contained_in(inner[4], outer[4]):
                        removed.add((inner[2], inner[3]))
            active.append(cur)

    kept = defaultdict(dict)
    for g, transcripts in transcript_dict.items():
        for t, c in transcripts.items():
            if (g, t) not in removed:
                kept[g][t] = c

    print("  Removed {} fully contained transcripts.".format(len(removed)))
    return kept


def _region_string(c):
    """UTR5/CDS/UTR3 as 1-based inclusive transcript coordinates, like the GENCODE fasta headers."""
    cds_s, cds_e = c["bed_CDS"]
    parts = []
    if cds_s > 0:
        parts.append("UTR5:1-{}".format(cds_s))
    parts.append("CDS:{}-{}".format(cds_s + 1, cds_e))
    if cds_e < c["length"]:
        parts.append("UTR3:{}-{}".format(cds_e + 1, c["length"]))
    return "|".join(parts)


def collapse_identical_transcripts(transcript_dict, groups_tsv):
    """
    Transcripts with identical merged exon structure (same chrom + strand) have
    identical RNA sequence, so reads cannot be assigned between them; in the
    reference they are different genes with different CDS on the same RNA
    (e.g. PYURF / PIGY). Keep only the member with the longest CDS, since RPFs
    come from the CDS. Every member is written to groups_tsv so the dropped ORF
    can still be quantified from coverage over the kept transcript.
    """
    by_structure = defaultdict(list)
    for g, transcripts in transcript_dict.items():
        for t, c in transcripts.items():
            key = (c["chrom"], c["strand"], tuple(map(tuple, _merge_intervals(c["exons"]))))
            by_structure[key].append((g, t, c))

    removed = set()
    n_groups = 0
    with open(groups_tsv, "w") as out:
        print("group_id\tmember_tx\tmember_gene\tregion\trepresentative", file=out)
        for members in by_structure.values():
            if len(members) < 2:
                continue
            n_groups += 1
            members.sort(key=lambda m: m[2]["bed_CDS"][1] - m[2]["bed_CDS"][0], reverse=True)
            group_id = "shared_" + "_".join(m[2].get("gene_name", m[0]) for m in members)
            for i, (g, t, c) in enumerate(members):
                print(group_id, t, c.get("gene_name", g), _region_string(c), i == 0,
                      sep="\t", file=out)
                if i > 0:
                    removed.add((g, t))

    kept = defaultdict(dict)
    for g, transcripts in transcript_dict.items():
        for t, c in transcripts.items():
            if (g, t) not in removed:
                kept[g][t] = c

    print("  {} identical-sequence groups; removed {} shorter-CDS members (see {}).".format(
        n_groups, len(removed), groups_tsv))
    return kept

#############################################################################################

def pick_protein_coding_transcripts( transcript_dict ):
    protein_coding_longest_transcripts = defaultdict(dict)

    for g, transcripts in transcript_dict.items():
        for t_name , t_contents in transcripts.items():
            if t_contents["gene_type"] == "protein_coding":
                protein_coding_longest_transcripts[g][t_name] = t_contents

    return protein_coding_longest_transcripts

#############################################################################################

#### CHANGE THIS FUNCTION AND USER IT

"""
def pick_transcripts_with_stop_codon( protein_coding_longest_transcripts, transcripts_fasta ):
    t_fasta_stream = FastaFile(transcripts_fasta)

    clipped = list()
    start_triplets = defaultdict(int)
    stop_triplets = defaultdict(int)
    empty_stop_triplets = defaultdict(dict)
    transcripts_with_stop_codon = defaultdict(dict)

    for this_fasta_entry in t_fasta_stream:
        contents = this_fasta_entry.header.strip().split("|")
        this_t = contents[0].split(".")[0]

        g_contents = contents[1].split(".")
        this_g = g_contents[0]

        if len(g_contents) >= 2 and "PAR_Y" in g_contents[1]:
            continue

        transcripts = protein_coding_longest_transcripts.get(this_g, None)
        if transcripts == None:
            continue
        if this_t in list(protein_coding_longest_transcripts[this_g].keys()):
            bed_CDS = protein_coding_longest_transcripts[this_g][this_t].\
                       get("bed_CDS", None)

            if bed_CDS != None:
                start_codon_coord = bed_CDS[0]
                start_triplets[ this_fasta_entry.sequence[ start_codon_coord:start_codon_coord+3 ]  ] += 1
                observed_stop_triplet = this_fasta_entry.sequence[ bed_CDS[1]:bed_CDS[1] + 3 ]
                stop_triplets[ observed_stop_triplet ] += 1
                if observed_stop_triplet == '':
                    empty_stop_triplets[this_g][this_t] = \
                       protein_coding_longest_transcripts[this_g][this_t]
                else:
                    transcripts_with_stop_codon[this_g][this_t] =\
                       protein_coding_longest_transcripts[this_g][this_t]

    return transcripts_with_stop_codon
"""

def pick_transcripts_with_proper_codons( appris_transcripts, fasta_file ):
    """
    Picks the transcripts with proper start and stop codon
    """
    proper_picks = defaultdict(dict)

    # Proper start codons are those that are 1 nucleotide away
    # from ATG (in Hamming Distance)

    proper_start_codons = ["ATG", "CTG", "GTG", "TTG",
                           "AAG", "ACG", "AGG",
                           "ATA", "ATC", "ATT"]

    proper_stop_codons = ['TAG', "TAA", "TGA"]

    t_fasta_stream = FastaFile(fasta_file)

    for this_fasta_entry in t_fasta_stream:
        contents = this_fasta_entry.header.strip().split("|")
        this_t   = contents[0].split(".")[0]

        g_contents = contents[1].split(".")
        this_g     = g_contents[0]

        transcripts = appris_transcripts.get(this_g, None)

        if transcripts == None:
            continue

        if this_t not in list(appris_transcripts[this_g].keys()):
            continue

        bed_CDS = appris_transcripts[this_g][this_t].get("bed_CDS", None)



        if bed_CDS != None:
            start_codon = this_fasta_entry.sequence[ bed_CDS[0] : bed_CDS[0] + 3 ]
            stop_codon  = this_fasta_entry.sequence[ bed_CDS[1] : bed_CDS[1] + 3 ]

            if not stop_codon:
                stop_codon  = this_fasta_entry.sequence[ bed_CDS[1] - 3 : bed_CDS[1] ]

            if start_codon in proper_start_codons and \
               stop_codon  in proper_stop_codons:
                proper_picks[this_g][this_t] = appris_transcripts[this_g][this_t]

    return proper_picks

#############################################################################################

def create_fasta_files( transcripts_fasta, transcripts_dict,
                        appris_selected_fasta, appris_discarded_fasta):

    t_fasta_stream = FastaFile(transcripts_fasta)
    selected_count = 0
    discarded_count = 0
    visited_transcripts = defaultdict(dict)

    with gzip.open(appris_selected_fasta, "wt") as selected_out_stream,\
         gzip.open(appris_discarded_fasta, "wt") as discarded_out_stream:

        for this_fasta_entry in t_fasta_stream:
            contents = this_fasta_entry.header.strip().split("|")
            this_t = contents[0].split(".")[0]
            g_contents = contents[1].split(".")
            this_g = g_contents[0]

            if len(g_contents) >= 2 and "PAR_Y" in g_contents[1]:
                continue

            if visited_transcripts.get(this_g, None) and \
               visited_transcripts[this_g].get(this_t, None):
               print("Warning: the following gene. "
                     "transcript pair appear more than once in the fasta file")
               print(">", this_g, this_t, "\n\n")

            visited_transcripts[this_g][this_t] = True

            if transcripts_dict.get( this_g, None ) and\
                transcripts_dict[this_g].get(this_t, None):

                print(this_fasta_entry, file=selected_out_stream)
                selected_count += 1

            else:
                print(this_fasta_entry, file=discarded_out_stream)
                discarded_count += 1

    return (selected_count, discarded_count)


#############################################################################################

def _prep_bed_entry( ref_name, start_end_pair, name ):
    pre_bed_list = [ ref_name, start_end_pair[0], start_end_pair[1],
                     name, 0, "+" ]
    bed_list = tuple(map( str, pre_bed_list ))
    return "\t".join( bed_list )

def _get_reg_region(ref_name, t_dict, region, bed_name ):
    result = None
    region_boundaries = t_dict.get(region, None)
    if region_boundaries != None:
        result = _prep_bed_entry( ref_name, region_boundaries, bed_name )
    return result


def _get_regions(ref_name, t_dict):
    entries = []

    for e in ( ("bed_UTR_5", UTR5_name), ("bed_CDS", CDS_name), ("bed_UTR_3", UTR3_name) ):
        this_entry = _get_reg_region(ref_name, t_dict, e[0], e[1] )
        if this_entry != None:
            entries.append(this_entry)

    return "\n".join(entries)


def _get_metagene_regions( ref_name, t_dict, metagene_radius ):
    result = ""

    bed_CDS = t_dict.get("bed_CDS", None)
    if bed_CDS == None:
        return ""

    # +1: Lst nucleotide is not included in the bed entry
    this_start_site = (bed_CDS[0] - metagene_radius , bed_CDS[0] + metagene_radius + 1 )
    this_stop_site = (bed_CDS[1] - metagene_radius, bed_CDS[1] + metagene_radius + 1 )

    result = _prep_bed_entry( ref_name, this_start_site, start_site_name )
    result += "\n" + _prep_bed_entry( ref_name, this_stop_site, stop_site_name )

    return result


def create_bed_files(transcripts_dict, selected_fasta, output_base):

    t_fasta_stream = FastaFile(selected_fasta)

    actual_regions_file = output_base + "_actual_regions.bed"

    myopen = open

    with myopen(actual_regions_file, "wt") as regions_out:

        for this_fasta_entry in t_fasta_stream:
            contents = this_fasta_entry.header.strip().split("|")
            this_t = contents[0].split(".")[0]

            g_contents = contents[1].split(".")
            this_g = g_contents[0]

            if len(g_contents) >= 2 and "PAR_Y" in g_contents[1]:
                continue

            transcripts = transcripts_dict.get(this_g, None)
            if transcripts == None:
                continue

            t_dict = transcripts.get(this_t, None)
            if t_dict == None:
                continue

            t_regions_bed = _get_regions(this_fasta_entry.header, t_dict)

            print( t_regions_bed, file = regions_out )


#############################################################################################

# Sanity Check Functions

def CDS_length( transcript_dict ):
    CDS_remainder_by_3 = [0,0,0]
    non_zero_remainder_transcripts = []

    for g, transcripts in transcript_dict.items():
        for t_name, t_contents in transcripts.items():
            bed_CDS = t_contents.get("bed_CDS", None)
            this_remainder = (bed_CDS[1] - bed_CDS[0]) % 3
            CDS_remainder_by_3[this_remainder] += 1
            if this_remainder != 0:
                non_zero_remainder_transcripts.append(t_name)

    return (CDS_remainder_by_3 , non_zero_remainder_transcripts)


"""
def start_stop_codon(transcripts_fasta, transcripts_dict):
    t_fasta_stream = FastaFile(transcripts_fasta)

    clipped = list()
    start_triplets = defaultdict(int)
    stop_triplets = defaultdict(int)
    empty_stop_triplets = defaultdict(dict)
    transcripts_with_stop_codon = defaultdict(dict)

    for this_fasta_entry in t_fasta_stream:
        contents = this_fasta_entry.header.strip().split("|")
        this_t = contents[0].split(".")[0]

        g_contents = contents[1].split(".")
        this_g = g_contents[0]
        if len(g_contents) >= 2 and "PAR_Y" in g_contents[1]:
            continue

        transcripts = transcripts_dict.get(this_g, None)

        if this_t in list(transcripts_dict[this_g].keys()):
            bed_CDS = transcripts_dict[this_g][this_t].get("bed_CDS", None)

            if bed_CDS != None:
                start_codon_coord = bed_CDS[0]
                start_triplets[ this_fasta_entry.sequence[ start_codon_coord:start_codon_coord+3 ]  ] += 1
                observed_stop_triplet = this_fasta_entry.sequence[ bed_CDS[1]:bed_CDS[1] + 3 ]
                stop_triplets[ observed_stop_triplet ] += 1
                if observed_stop_triplet == '':
                    empty_stop_triplets[this_g][this_t] = transcripts_dict[this_g][this_t]
                else:
                    transcripts_with_stop_codon[this_g][this_t] = transcripts_dict[this_g][this_t]

    start_sorted = sorted(start_triplets.items(), key=lambda kv: kv[1])
    start_sorted.reverse()

    stop_sorted = sorted(stop_triplets.items(), key=lambda kv: kv[1])
    stop_sorted.reverse()

    return (start_sorted, stop_sorted)
"""

def find_cds_lengths(transcripts_dict):
    cds_lengths = []
    t_lengths = []

    for g, transcripts in transcripts_dict.items():
        for t_name, t_contents in transcripts.items():
            cds_lengths.append( t_contents["bed_CDS"][1] - t_contents["bed_CDS"][0] )
            t_lengths.append(t_contents["length"])

    return (cds_lengths, t_lengths)

def create_transcript_lengths_file(input_fasta_file, output_tsv_file):
    with FastaFile(input_fasta_file) as fasta_stream,\
         open(output_tsv_file, "w") as output_stream:

         for fasta_entry in fasta_stream:
             this_line = "{}\t{}".format( fasta_entry.header,
                                          len(fasta_entry.sequence)  )
             print(this_line, file = output_stream)

#############################################################################################

# Note we are using gene and transcript interchangably. We are basically picking
# one transcript from each gene and filter this list.

def main():

    arguments = get_arguments()
    gtf_file = arguments.g
    transcript_fasta_file = arguments.f
    appris_transcript_file = arguments.a
    output_base = arguments.o

    print("Reading the gtf file {}".format(gtf_file))
    gtf_all = get_gtf_contents_genome( gtf_file )
    find_rel_utr_cds_positions(gtf_all)
    print("Done reading gtf file.")

    appris_principal_transcripts, mane_select_transcripts = \
       read_principal_appris_trascript_list(appris_transcript_file,
                                            gtf_all)
    print("Started with {} principal appris transcripts".format( len(appris_principal_transcripts) ))

    # Step 1: Pick MANE_Select+PRINCIPAL transcript, or longest PRINCIPAL as fallback
    longest_appris_transcripts = pick_best_appris_transcripts(
        appris_principal_transcripts, mane_select_transcripts)
    print("Picked best transcript for each gene. Total number of picks: ",
           len(longest_appris_transcripts) )

    # Step 2: Pick only protein coding genes
    protein_coding_transcripts = pick_protein_coding_transcripts(longest_appris_transcripts)
    print("Picked protein coding transcripts. Total number of picks: ",
           len(protein_coding_transcripts))

    # Step 3: Remove read-through genes (tagged "readthrough_gene" in the GTF)
    readthrough_gene_ids = read_readthrough_gene_ids(gtf_file)
    print("Found {} read-through genes in the gtf file.".format(len(readthrough_gene_ids)))
    non_readthrough_transcripts = \
        remove_readthrough_genes(protein_coding_transcripts, readthrough_gene_ids)
    print("Removed read-through genes. Total number of picks: ",
           len(non_readthrough_transcripts))

    # Step 4: Pick transcripts with stop codon
    transcripts_with_proper_codons = \
        pick_transcripts_with_proper_codons(non_readthrough_transcripts, transcript_fasta_file)
    print("Picked transcripts with stop codon. Total number of picks: ",
           len(transcripts_with_proper_codons))

    # Step 5: Remove transcripts fully contained in a longer selected transcript
    transcripts_with_proper_codons = \
        remove_contained_transcripts(transcripts_with_proper_codons)
    print("Removed fully contained transcripts. Total number of picks: ",
           len(transcripts_with_proper_codons))

    # Step 6: Collapse transcripts with identical sequence; keep the longest CDS,
    # record all members in the sequence-groups side table
    transcripts_with_proper_codons = \
        collapse_identical_transcripts(transcripts_with_proper_codons,
                                       output_base + "_sequence_groups.tsv")
    print("Collapsed identical-sequence transcripts. Total number of picks: ",
           len(transcripts_with_proper_codons))

    # Create Fasta Files
    appris_selected_fasta = output_base + "_selected.fa.gz"
    appris_discarded_fasta = output_base + "_discarded.fa.gz"

    print("Writing fasta files.")

    selected_count, discarded_count = \
       create_fasta_files( transcript_fasta_file, transcripts_with_proper_codons,
                           appris_selected_fasta, appris_discarded_fasta )
    print("Fasta Files written:\n{} transcripts selected.".format(selected_count))
    print("{} transcripts discarded.".format(discarded_count))

    # Create Bed Files

    print("Creating bed files.")
    create_bed_files(transcripts_with_proper_codons, appris_selected_fasta, output_base)

    print("Creating transcript lengths file.")
    transcript_length_file = output_base + "_transcript_lengths.tsv"
    create_transcript_lengths_file(appris_selected_fasta, transcript_length_file)


    ## Display sanity checks
    ## CDS Length
    remainder_array, transcripts_with_CDS_violation = CDS_length(transcripts_with_proper_codons)
    print("There are {} transcripts with CDS length divisible by 3".format(remainder_array[0]))
    CDS_violators_count = remainder_array[1] + remainder_array[2]
    print( "There are", CDS_violators_count, "transcripts with CDS length NOT divisible by 3." )
    if CDS_violators_count > 0:
        print("Some CDS violators are:")
        print(transcripts_with_CDS_violation[:min(10, CDS_violators_count)])


    cds_lengths, t_lengths = find_cds_lengths(transcripts_with_proper_codons)
    print("Median CDS length is ", statistics.median(cds_lengths))
    print("Median transcript length is ", statistics.median(t_lengths))


if __name__ == "__main__":
    main()
