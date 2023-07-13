#!python3

import argparse
from pathlib import Path
import sys
import math
import random
import statistics

VERSION = "1.0.0"

def main():

    # Get the options passed to the program
    options = read_options()

    # Read in the codon table for the species they're using. We get back
    # the translation table to be able to make protein sequence, the full
    # frequency values for codon usage, and the w values which are the 
    # frequencies normalised to the most commonly used codon
    codon_table, amino_acid_frequencies, w_values = load_codon_table(options)


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
    for seq_name,cds_sequence in read_cds_sequence(options.seqfile):
        log("Processing "+seq_name)

        if options.gc is None:
            gc = calculate_percent_gc(cds_sequence)
        else:
            gc = options.gc

        weighted_codons = calculate_weighted_codons(amino_acid_frequencies, gc)

        # We translate the CDS to get the protein sequence.
        try:
            protein_sequence = translate_cds(cds_sequence, codon_table)
        except:
            continue
        

        # We calculate the true observed CAI for this sequence
        log("Calcuating true CAI")
        true_cai = calculate_cai(cds_sequence, protein_sequence, w_values)

        # We generate a set of random sequences based on the composition of the
        # real protein and using the GC content of the organism to pick which 
        # codon to use from the set of synonymous options.  We then calculate the
        # CAI values from these to generate a background
        log("Generating background CAI distribution")
        background_cai = generate_background_cai(protein_sequence, weighted_codons, w_values, options)

        # We now judge the true CAI against the set of random sequences.
        log("Calculating CAI statistics")
        expected_cai = statistics.mean(background_cai)
        cai_stdev = statistics.stdev(background_cai)
        normalised_cai = true_cai/expected_cai
        cai_zscore = (true_cai-expected_cai)/cai_stdev

        results = [
            seq_name,
            len(protein_sequence),
            gc,
            true_cai,
            expected_cai,
            cai_stdev,
            normalised_cai,
            cai_zscore,
            ":".join([str(x) for x in background_cai])
        ]

        print("\t".join([str(x) for x in results]), file=out)


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

def generate_random_sequence(input_seq, weighted_codons):
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

    for i in range(options.samples):
        debug("Generating random sequence "+str(i+1))
        seq = generate_random_sequence(protein_sequence, weighted_codons)
        background_cai.append(calculate_cai(seq,protein_sequence,w_values))

    return background_cai

def calculate_cai(cds, aa, w_values):

    # Calculation comes from https://journals.sagepub.com/doi/10.1177/117693430700300028

    # First we count the observed codons
    codon_count = {}
    for i in range(len(aa)):
        codon = cds[i:i+3]
        if not codon in codon_count:
            codon_count[codon] = 0
        
        codon_count[codon] += 1

    # Now we get the sum of the ln(w) values for each codon
    cai_sum = 0
    for codon,count in codon_count.items():
        cai_sum += count * math.log(w_values[codon])

    debug("Log w sum is"+str(cai_sum))
    # We divide the sum by the number of amino acids to 
    # get the mean value per amino acid

    cai_sum /= len(aa)

    debug("Mean log w is"+str(cai_sum))

    # Finally we get the cai value by taking the exponent of
    # the average value

    cai = math.e ** cai_sum
    debug("CAI is "+str(cai))

    return cai


def translate_cds(cds, codons):
    log("Translating CDS to AA sequence")
    aa_seq = ""

    index = 0
    while(index+1 < len(cds)):
        aa_seq += codons[cds[index:index+3]]

        index += 3

    if aa_seq.endswith("*"):
        aa_seq = aa_seq[:-1]

    if "*" in aa_seq:
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

    w_values = {}

    for data in amino_acids.values():
        data.sort(key=lambda x:x["freq"], reverse=True)
        highest = data[0]["freq"]
        for item in data:
            w_values[item["codon"]] = item["freq"]/highest

    return codons, amino_acids, w_values


def read_cds_sequence(seqfile):
    log("Reading sequences from"+str(seqfile))

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
            return (seq_name,sequence)

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
    parser.add_argument("--gc", default=None, help="Manually specific GC content (uses sequence GC otherwise)")
    parser.add_argument("--quiet", action="store_true", help="Suppress all progress messages")
    parser.add_argument("--debug", action="store_true", help="Show verbose debugging messages")
    parser.add_argument("--random", type=str, default="markov", help="Method to generate random sequences, values are 'markvov' (default) or 'shuffle'")

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