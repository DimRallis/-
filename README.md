# KAMY
KAMY is an pipeline that uses the kmer-counting algorithm Meryl to identify sex-specific regions in an assembly. Users should install an older version of Meryl (https://kmer.sourceforge.net) and add the KAMY.sh, split_genome_to_Ns.pl, fasta2fake_fastq.pl and kmer_match_genome_sizes.pl scripts to their $PATH.

The KAMY.sh runs as an automated script. It takes as input seperate male and female paired-end sequencing reads and a set of query nucleotide sequences (chromosomes, scaffolds, transcripts, etc). Users are welcome to modify the KAMY.sh script according to their needs. Sequencing files should be in .fastq format and nucleotide sequence file should include a .fasta extension.

Once KAMY is run over a set of paired-end sequencing files, the respective kmer databases are saved in the output directory and the query set is assessed for Male/Female/Common kmer coverage. After that, users can save time and computational resources when running KAMY against other query datasets, by calling the -skip_counting function. This offers the advantage of rapidly determining male-specific regions during genomic projects, since different assembly versions or transcriptomic datasets can be instantly assessed for Y-linkage without the need for de novo kmer-counting.

KAMY should be run as:

bash KAMY.sh [OPTIONS]

Use -h to see the full list of options.
