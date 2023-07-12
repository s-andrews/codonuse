#!python3

import argparse

def main():
    options = read_options()
    codon_table = load_codon_table(options)
    background_cai = generate_background_cai(options)
    true_cai = calculate_cai(options)

    write_results()


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


if __name__ == "__main__":
    main()