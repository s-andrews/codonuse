# codonuse - Efficient calculation of codon optimisation

## Background
Although the genetic code contains redundancy, with the same amino acid often being able to be encoded by multiple codons, the distribution of tRNAs in an organism means that some codons are used more efficiently than others, and that by changing from an optimised to a non-optimised codon balance you can potentially modulate the rate at which a transcript is traslated.

It is useful to have a method to quantitate the degree to which the codon usage for a particular sequence is optimised and that is what this program does.  It uses the method described in https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-65 to produce a number of metrics to assess the optimisation of codon usage.

The basic measurement used is the Codon Adaption Index, which was originall described in https://academic.oup.com/nar/article/15/3/1281/1166844, but with a more efficient implemenation described in https://journals.sagepub.com/doi/10.1177/117693430700300028.  This program calculates the CAI for a sequence of interest, but then scales this against a series of pseudo-random sequences which randomly select from the available codons based solely on the GC content of their third base, so that biases in genome composition don't affect the results.

## Running the program

The structure for running the program is as follows

```
usage: codonuse.py [-h] [--version] [--outfile OUTFILE] [--samples SAMPLES] [--gc GC] [--quiet] [--debug]
                   [--random {markov,shuffle}] [--threads THREADS]
                   species seqfile

Calculate the Codon Adaptive Index for a sequence and compare to a reference set

positional arguments:
  species               Name of species - must match a codon file in the 'tables' directory
  seqfile               Filename for (multi-)fasta format file of mRNA coding sequence

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --outfile OUTFILE     Name of file into which to write results (default codonuse_output.txt)
  --samples SAMPLES     Number of random sequences to generate (default 500)
  --gc GC               Manually specific GC content (uses sequence GC otherwise)
  --quiet               Suppress all progress messages
  --debug               Show verbose debugging messages
  --random {markov,shuffle}
                        Method to generate random sequences, values are 'markvov' (default) or 'shuffle'
  --threads THREADS     Number of threads to use

Report problems at https://github.com/s-andrews/codonuse/issues
```

You need to supply 

1. A fasta format sequence file of coding sequences.  This can contain multiple entries and they will all be analysed.  The sequences must only be the CDS sequences (with or without stop codon is OK), with no 5' or 3' UTR.

2. A species name - this must match the name of one of the codon preference tables in the ```tables``` directory of the program.  We have included some common species but you can add support for more by adding in your own codon table based on the data in https://www.kazusa.or.jp/codon/

To run the program with a single example sequence you can run:

```
python3 codonuse.py drosophila test_data/wingless.fa
```

The output will be written to a file called ```codonuse_output.txt``` 

For a more complete analysis you can run:

```
python3 codonuse.py --outfile test_data/drosophila_cai.txt --gc=53 --threads=4 drosophila test_data/drosophila.fa
```

This will generate output for a larger collection of protein coding genes in the drosophila genome.



## Output Format

The output is a tab-delimited text file with one line of results for each sequence analysed.  The columns in the output are:

1. Sequence Name
2. Protein Length
3. GC (either the %GC of the sequence, or the value provided in the ```--gc``` argument)
4. The CAI value for the sequence
5. The mean CAI value for the random sequences
6. The stdev of CAI values for the random sequences
7. The normalised CAI (true/mean random)
8. The normalised z-score ((true-mean)/stdev)
9. A colon delimited list of the CAI scores for the background set
