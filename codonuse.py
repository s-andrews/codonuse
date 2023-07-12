#!python3

import argparse
from pathlib import Path
import sys
import math

def main():
    options = read_options()

    seq_name,cds_sequence = read_cds_sequence(options)
    codon_table, amino_acid_frequencies, w_values = load_codon_table(options)
    protein_sequence = translate_cds(cds_sequence, codon_table)
    true_cai = calculate_cai(cds_sequence, protein_sequence, w_values)
    background_cai = generate_background_cai(options)

    write_results()


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


def read_cds_sequence(options):
    log("Reading CDS from"+options.seqfile)
    with open(options.seqfile,"rt", encoding="utf8") as infh:
        header = infh.readline()
        if not header.startswith(">"):
            raise Exception("Sequence file"+options.seqfile+" isn't in fasta format (first line didn't start with >)")

        seqname = header.split()[0][1:]

        debug("Found sequence '"+seqname+"'")

        sequence = ""
        for line in infh:
            line = line.strip()
            line = line.replace(" ","")
            line = line.upper()
            line = line.replace("U","T")
            sequence += line

        if not len(sequence)%3 == 0:
            raise Exception("Length of sequence ("+str(len(sequence))+") was not a multiple of 3 - doesn't look like coding sequence")

        debug("Read sequence: "+sequence)

        return(seqname, sequence)


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

    parser.add_argument("seqfile", type=str, help="Filename for fasta format file of mRNA coding sequence")
    parser.add_argument("species", type=str, help="Name of species - must match a codon file in the 'tables' directory")
    parser.add_argument("--samples", type=int, default=500, help="Number of random sequences to generate")
    parser.add_argument("--gc", default=None, help="Manually specific GC content (uses sequence GC otherwise)")
    parser.add_argument("--quiet", action="store_true", help="Suppress all progress messages")
    parser.add_argument("--debug", action="store_true", help="Show verbose debugging messages")

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