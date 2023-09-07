#!/usr/bin/env python3
#
# This program is free software: you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free 
# Software Foundation, either version 3 of the License, or (at your option) 
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along 
# with this program. If not, see <https://www.gnu.org/licenses/>. 

import argparse
from pathlib import Path
import sys
import math
import random
import statistics
import hashlib
from multiprocess import Pool

VERSION = "1.1.0"

QUIET = False
DEBUG = False

def main():

    # Get the options passed to the program
    options = read_options()

    # Read in the codon table for the species they're using. We get back
    # the translation table to be able to make protein sequence, the full
    # frequency values for codon usage, and the w values which are the 
    # frequencies normalised to the most commonly used codon
    codon_table, amino_acid_frequencies, w_values = load_codon_table(options)

    # Make a weighted set of codons.  We can only do this if they specified
    # a fixed GC for all analyses, otherwise we have to recalculate per
    # sequence as the GC will be different in all cases.
    weighted_codons = None

    if options.gc is not None:
        weighted_codons = calculate_weighted_codons(amino_acid_frequencies, options.gc)

    # Open the output file and print the headers.
    result_names = [
        "seq_name",
        "protein_length",
        "gc",
        "true_cai",
        "background_mean_cai",
        "background_cai_stdev",
        "normalised_cai",
        "cai_zscore",
        "background_values"
    ]

    out = open(options.outfile, "wt", encoding="utf8")

    print("\t".join(result_names), file=out)

    # The seqfile can contain multiple sequences so we'll iterate through them
    # We read in the input sequence and get back both the name and the DNA
    # sequence.  Sequence will be DNA (not RNA) and upper case

    with Pool(options.threads) as pool:
        results = pool.map(lambda x: process_sequence(x[0],x[1],options,amino_acid_frequencies,codon_table,w_values,weighted_codons),read_cds_sequence(options.seqfile))


    for r in results:
        if r is not None: # We'd get this if the translation failed
            print(r,file=out)

    # for seq_name,cds_sequence in read_cds_sequence(options.seqfile):
    #     process_sequence(seq_name,cds_sequence,options,amino_acid_frequencies,codon_table,w_values,out)


def process_sequence(seq_name, cds_sequence,options,amino_acid_frequencies,codon_table,w_values,weighted_codons):
    log("Processing "+seq_name)

    if options.gc is None:
        gc = calculate_percent_gc(cds_sequence)
        weighted_codons = calculate_weighted_codons(amino_acid_frequencies, gc)
    else:
        gc = options.gc

    
    # We translate the CDS to get the protein sequence.
    try:
        protein_sequence = translate_cds(cds_sequence, codon_table)
    except:
        return
    

    # We calculate the true observed CAI for this sequence
    debug("Calcuating true CAI")
    true_cai = calculate_cai(cds_sequence, protein_sequence, w_values)

    # We generate a set of random sequences based on the composition of the
    # real protein and using the GC content of the organism to pick which 
    # codon to use from the set of synonymous options.  We then calculate the
    # CAI values from these to generate a background
    debug("Generating background CAI distribution")

    background_cai = generate_background_cai(protein_sequence, weighted_codons, w_values, options)

    # We now judge the true CAI against the set of random sequences.
    log("Calculating CAI statistics")
    expected_cai = statistics.mean(background_cai)
    cai_stdev = statistics.stdev(background_cai)

    # In corner cases it's possible for the the CAI values in the background
    # to all be the same, generating a zero stdev, so in these cases we have
    # to set an arbitrarily low value so it doesn't generate a divide by zero
    # error
    if cai_stdev == 0:
        cai_stdev = 0.001
        log("Setting zero stdev to 0.001")

    normalised_cai = true_cai/expected_cai
    cai_zscore = (true_cai-expected_cai)/cai_stdev

    results = [
        seq_name,
        len(protein_sequence),
        gc,
        round(true_cai,3),
        round(expected_cai,3),
        round(cai_stdev,3),
        round(normalised_cai,3),
        round(cai_zscore,3),
        ":".join([str(round(x,3)) for x in background_cai])
    ]

    return "\t".join([str(x) for x in results])


def calculate_percent_gc(seq):
    log("Calculating GC from sequence as manual value not specified")
    gc = 0

    for i in seq:
        if i=="G" or i=="C":
            gc += 1

    gc /= len(seq)
    gc *= 100

    return gc

def calculate_weighted_codons (aafreqs, gc_percent):

    log("Calculating weighted codons")

    # The weighting values are just based on the GC content of the last base.

    weighted_codons = {}

    for aa in aafreqs.keys():
        weighted_codons[aa] = {"codons": [], "weights": []}
        for i in aafreqs[aa]:
            codon = i["codon"]
            weighted_codons[aa]["codons"].append(codon)
            if codon.endswith("G") or codon.endswith("C"):
                weighted_codons[aa]["weights"].append(gc_percent)
            else:
                weighted_codons[aa]["weights"].append(100-gc_percent)


    return weighted_codons

def backtranslate(seq, codons):
    dna = []

    for aa in seq:
        dna.append(random.choices(codons[aa]["codons"], weights=codons[aa]["weights"])[0])

    return "".join(dna)

def generate_random_sequence(input_seq, weighted_codons, method):
    # If we're not doing markov sampling then we just reuse the existing protin
    # sequence and reassign codons
    protein_seq = input_seq

    if method=="markov":
        protein_seq = []

        for _ in range(len(input_seq)):
            protein_seq.append(random.choice(input_seq))

        protein_seq = "".join(protein_seq)
    debug("Random protein "+protein_seq)

    dna_seq = backtranslate(protein_seq, weighted_codons)
    debug("Random DNA "+dna_seq)
    return dna_seq

def generate_background_cai(protein_sequence, weighted_codons, w_values, options):
    background_cai = []

    # To make this deterministic we set the random seed
    # to a value based on the protein sequence so we'll
    # always get the same answer for the same sequence

    random.seed(
        int(
            hashlib.md5(protein_sequence.encode("utf8")).hexdigest(),
            base=16
        )
    )

    for i in range(options.samples):
         debug("Generating random sequence "+str(i+1))
         background_cai.append(generate_random_cai(protein_sequence,weighted_codons,options.random,w_values))

    return background_cai


def generate_random_cai(protein_sequence, weighted_codons, randomtype, w_values):
    seq = generate_random_sequence(protein_sequence, weighted_codons, randomtype)
    cai = calculate_cai(seq,protein_sequence,w_values)
    return cai

def calculate_cai(cds, aa, w_values):

    # Calculation comes from https://journals.sagepub.com/doi/10.1177/117693430700300028

    # We need to get the frequency of codon uses for each amino acid
    aa_count = {}
    for i in range(len(aa)):
        thisaa = aa[i]
        codon = cds[i*3:(i*3)+3]
        if not thisaa in aa_count:
            aa_count[thisaa] = {}

        if not codon in aa_count[thisaa]:
            aa_count[thisaa][codon] = 0

        aa_count[thisaa][codon] += 1

    # We need to convert the counts into frequencies per amino acid
    total = 0
    for aa in aa_count.keys():
        for count in aa_count[aa].values():
            total += count

    for aa in aa_count.keys():
        for codon in aa_count[aa].keys():
            aa_count[aa][codon] /= total


    total_freqs = 0
    total_logw = 0

    breakpoint()
    for aa in aa_count.keys():
        for codon in aa_count[aa].keys():
            if not codon in w_values:
                # This must come from a single codon family so skip it
                continue
            total_freqs += aa_count[aa][codon]
            total_logw += math.log(w_values[codon])*aa_count[aa][codon]

    debug(f"Total logw {total_logw} total freqs {total_freqs}")

    cai = math.e ** (total_logw/total_freqs)

    log(f"CAI is {cai}")
    return cai


def translate_cds(cds, codons):
    debug("Translating CDS to AA sequence")
    aa_seq = ""

    index = 0
    while(index+1 < len(cds)):
        aa_seq += codons[cds[index:index+3]]

        index += 3

    if aa_seq.endswith("*"):
        aa_seq = aa_seq[:-1]

    if "*" in aa_seq:
        log("Sequence translated with internal stops. Can't continue")
        raise Exception("Sequence translated with internal stops. Can't continue")


    debug("Translated protein: "+aa_seq)

    return aa_seq


def load_codon_table(options):
    species = options.species
    log("Reading coding table for "+species)
    filename = species+".txt"
    file = Path(__file__).parent / "tables" / filename

    debug("Codon file is "+str(file))

    codons = {}
    amino_acids = {}

    with open(file) as infh:
        for line in infh:
            line = line.strip()

            # We need to remove spaces from things like "(  1234)" otherwise it breaks
            # the splitting on whitespace
            while "( " in line:
                line = line.replace("( ","(")

            sections = line.split()
            start = 0

            # The lines consist of multiple sections which come in groups of 5. We only
            # need the first 3 (codon, aa, freq) from each set of 5
            while (start+1 < len(sections)):
                codon = sections[start]
                codon = codon.upper().replace("U","T")
                aa = sections[start+1].upper()
                freq = float(sections[start+2])

                debug("Found codon "+codon+" for aa "+aa+" with freq "+str(freq))

                codons[codon] = aa

                if not aa in amino_acids:
                    amino_acids[aa] = []

                amino_acids[aa].append({"codon":codon, "freq":freq})

                start += 5

    # We need to convert the frequency values into a w value, which is
    # the frequency expressed as a proportion of the highest frequency
    # so the most common codon gets a value of 1 and the others go down
    # from there.

    # A complication is that the variation is only in the last base of 
    # the triplet and that some amino acids have more than one 'family'
    # where there is more than one combination of first codons which 
    # encode the same amino acid.  We need to calculate the w values
    # separately for each of these families.

    w_values = {}    

    for aa,data in amino_acids.items():

        # First split into families
        families = {}
        for item in data:
            family = item["codon"][:2]
            if not family in families:
                families[family] = []
            families[family].append(item)

        # Now calculate w values for each family
        for family in families.values():
            # If we have a family with only one member then we exclude it as it 
            # provides no useful information
            if len(family) == 1:
                debug(f"Skipping family {family[0]['codon'][:2]} with only 1 member")
                continue
            family.sort(key=lambda x:x["freq"], reverse=True)
            highest = family[0]["freq"]
            debug(f"For AA {aa} family {family[0]['codon'][:2]} highest frequency is {highest}")
            for item in family:
                w_values[item["codon"]] = item["freq"]/highest
                debug(f"For codon {item['codon']} freq is {item['freq']} w is {w_values[item['codon']]}")

    return codons, amino_acids, w_values


def read_cds_sequence(seqfile):
    log("Reading sequences from "+str(seqfile))

    seq_name = ""
    sequence = ""

    with open(seqfile,"rt", encoding="utf8") as infh:

        for line in infh:
            if line.startswith(">"):
                if seq_name:
                    yield(seq_name,sequence)

                seq_name = line.split()[0][1:]
                sequence = ""

            else:
                line = line.strip()
                line = line.replace(" ","")
                line = line.upper()
                line = line.replace("U","T")
                sequence += line

        if seq_name:
            yield (seq_name,sequence)

def log(message):
    if not QUIET:
        print("LOG:",message, file=sys.stderr)

def debug(message):
    if DEBUG:
        print("DEBUG:",message, file=sys.stderr)


def read_options():

    parser = argparse.ArgumentParser(
        prog="codonuse.py",
        description="Calculate the Codon Adaptive Index for a sequence and compare to a reference set",
        epilog="Report problems at https://github.com/s-andrews/codonuse/issues"
    )

    parser.add_argument('--version', action='version', version=f"codonuse.py version {VERSION}")
    parser.add_argument("species", type=str, help="Name of species - must match a codon file in the 'tables' directory")
    parser.add_argument("seqfile", type=str, help="Filename for (multi-)fasta format file of mRNA coding sequence")
    parser.add_argument("--outfile", type=str, help="Name of file into which to write results (default codonuse_output.txt)", default="codonuse_output.txt")
    parser.add_argument("--samples", type=int, default=500, help="Number of random sequences to generate (default 500)")
    parser.add_argument("--gc", type=int, default=None, help="Manually specific GC content (uses sequence GC otherwise)")
    parser.add_argument("--quiet", action="store_true", help="Suppress all progress messages")
    parser.add_argument("--debug", action="store_true", help="Show verbose debugging messages")
    parser.add_argument("--random", type=str, choices=['markov','shuffle'], default="markov", help="Method to generate random sequences, values are 'markvov' (default) or 'shuffle'")
    parser.add_argument("--threads", type=int, help="Number of threads to use", default=1)

    options = parser.parse_args()

    global QUIET
    QUIET = options.quiet

    global DEBUG
    DEBUG = options.debug
    if DEBUG:
        QUIET = False


    return(options)


if __name__ == "__main__":
    main()
