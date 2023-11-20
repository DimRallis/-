#!/bin/bash

#Save the split_genome_to_Ns.pl fasta2fake_fastq.pl and kmer_match_genome_sizes.pl scrips to your bin file

print_help() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -skip_counting           Skip kmer counting and database generation if previously done, users should specify after this option the directory under which the kmer database lives (optional)"
    echo "  -skip_splitting          Do not split the input genome into contigs (optional)"
    echo "  -threads                 Number of threads, recommended to use >20 (default: 1)"
    echo "  -male_1                  Path to male read 1 fastq file"
    echo "  -male_2                  Path to male DNA read 2 fastq file"
    echo "  -female_1                Path to female DNA read 1 fastq file"
    echo "  -female_2                Path to female DNA read 2 fastq file"
    echo "  -genome                  Path to the genome fasta file"
    echo "  -out_dir                 Output directory"
    echo "  -h                       Print this help message"
    exit 0
}

SKIP_COUNTING=false
DB_DIR=""
SKIP_SPLITTING=false
THREADS=1


while [[ $# -gt 0 ]]; do
    case "$1" in
        -skip_counting)
            SKIP_COUNTING=true
            shift
            DB_DIR="$1"
            shift
            ;;
        -skip_splitting)
            SKIP_SPLITTING=true
            shift
            ;;
        -threads)
            shift
            THREADS=$1
            shift
            ;;
        -male_1)
            shift
            MALE_READ_1=$1
            shift
            ;;
        -male_2)
            shift
            MALE_READ_2=$1
            shift
            ;;
        -female_1)
            shift
            FEMALE_READ_1=$1
            shift
            ;;
        -female_2)
            shift
            FEMALE_READ_2=$1
            shift
            ;;
        -genome)
            shift
            GENOME_INPUT=$1
            shift
            ;;
        -out_dir)
            shift
            OUTPUT_DIR=$1
            shift
            ;;
        -h)
            print_help
            ;;
        *)
            echo "Error: Unknown option or incorrect number of arguments."
            print_help
            ;;
    esac
done

if [ "$SKIP_COUNTING" = true ] && [ -f "$DB_DIR/Male.kmer21.cut2.mcdat" ] && [ -f "$DB_DIR/Female.kmer21.cut2.mcdat" ] && [ -f "$DB_DIR/Common.kmer21.cut2.mcdat" ]; then
    echo "mcdat files exist"

elif [ ! -f "$MALE_READ_1" ] || [ ! -f "$MALE_READ_2" ] || [ ! -f "$FEMALE_READ_1" ] || [ ! -f "$FEMALE_READ_2" ] || [ ! -f "$GENOME_INPUT" ]; then
    echo "Error: One or more input files do not exist."
    print_help
fi

if [ "$SKIP_COUNTING" = true ] && [ -f "$GENOME_INPUT" ] && [ ! -z "$OUTPUT_DIR" ]; then
    echo "No read files needed, reading from existing kmer databases"


    if [ ! -d "$OUTPUT_DIR" ]; then
      mkdir -p "$OUTPUT_DIR" || { echo "Error: Unable to create output directory $OUTPUT_DIR"; exit 1; }
    fi

    if [ "$SKIP_SPLITTING" = false ]; then

      echo "Splitting genome to Ns"
      split_genome_to_Ns.pl $GENOME_INPUT
      GENOME="${GENOME_INPUT}.split"
      GENOME_INDEX="${GENOME}.fai"
      FAKE_FASTQ=$(echo $GENOME | sed 's/\.split$/.fastq/')

      if [ "$OUTPUT_DIR" != ./ ]; then
      mv "$GENOME" "$OUTPUT_DIR/$GENOME"
      fi

    else
      echo "Genome will not be split"
      GENOME=$GENOME_INPUT
      GENOME_INDEX="${GENOME}.fai"
      FAKE_FASTQ=$(echo $GENOME | sed 's/\.fasta$/.fastq/')
    fi

    echo "Generating fake fastq genome file"
    fasta2fake_fastq.pl $GENOME

    if [ "$OUTPUT_DIR" != ./ ]; then
    mv "$FAKE_FASTQ" "$OUTPUT_DIR/$FAKE_FASTQ"
    fi

    echo "Calculating max sequence length"
    MAX_LENGTH=$(awk -F'\t' -v col="2" '{ if ($col > max) max = $col } END { print max }' "Genome_sizes.tsv")
    ROUNDED_MAX_LENGTH=$(echo "scale=0; ((($MAX_LENGTH + 99) / 100) * 100)" | bc)

    echo "Indexing the genome"
    samtools faidx $GENOME

    if [ "$OUTPUT_DIR" != ./ ]; then
    mv $GENOME_INDEX $OUTPUT_DIR/$GENOME_INDEX
    fi

    cut -f1,2 $OUTPUT_DIR/$GENOME_INDEX > $OUTPUT_DIR/Genome_sizes.tsv

    if [ "$SKIP_COUNTING" = false ]; then

      echo "Splitting input read files into smaller chunks"
      SEQ_COUNT_MALE_1=$(echo $(cat $MALE_READ_1|wc -l)/4|bc)
      SEQ_PER_FILE_MALE_1=$(echo "($SEQ_COUNT_MALE_1 / 2) + 1" | bc)
      LINES_PER_FILE_MALE_1=$(echo "$SEQ_PER_FILE_MALE_1 * 4" | bc)

      SEQ_COUNT_MALE_2=$(echo $(cat $MALE_READ_2|wc -l)/4|bc)
      SEQ_PER_FILE_MALE_2=$(echo "($SEQ_COUNT_MALE_2 / 2) + 1" | bc)
      LINES_PER_FILE_MALE_2=$(echo "$SEQ_PER_FILE_MALE_2 * 4" | bc)

      SEQ_COUNT_FEMALE_1=$(echo $(cat $FEMALE_READ_1|wc -l)/4|bc)
      SEQ_PER_FILE_FEMALE_1=$(echo "($SEQ_COUNT_FEMALE_1 / 2) + 1" | bc)
      LINES_PER_FILE_FEMALE_1=$(echo "$SEQ_PER_FILE_FEMALE_1 * 4" | bc)

      SEQ_COUNT_FEMALE_2=$(echo $(cat $FEMALE_READ_2|wc -l)/4|bc)
      SEQ_PER_FILE_FEMALE_2=$(echo "($SEQ_COUNT_FEMALE_2 / 2) + 1" | bc)
      LINES_PER_FILE_FEMALE_2=$(echo "$SEQ_PER_FILE_FEMALE_2 * 4" | bc)


      split -l $LINES_PER_FILE_MALE_1 $MALE_READ_1 $OUTPUT_DIR/${MALE_READ_1}_
      split -l $LINES_PER_FILE_MALE_2 $MALE_READ_2 $OUTPUT_DIR/${MALE_READ_2}_
      split -l $LINES_PER_FILE_FEMALE_1 $FEMALE_READ_1 $OUTPUT_DIR/${FEMALE_READ_1}_
      split -l $LINES_PER_FILE_FEMALE_2 $FEMALE_READ_2 $OUTPUT_DIR/${FEMALE_READ_2}_

      echo "Generating kmer table from reads"
      meryl -B -C -m 21 -s $OUTPUT_DIR/${MALE_READ_1}_aa -o $OUTPUT_DIR/Male_R1_aa.kmer21 -threads $THREADS
      meryl -B -C -m 21 -s $OUTPUT_DIR/${MALE_READ_1}_ab -o $OUTPUT_DIR/Male_R1_ab.kmer21 -threads $THREADS

      meryl -B -C -m 21 -s $OUTPUT_DIR/${MALE_READ_2}_aa -o $OUTPUT_DIR/Male_R2_aa.kmer21 -threads $THREADS
      meryl -B -C -m 21 -s $OUTPUT_DIR/${MALE_READ_2}_ab -o $OUTPUT_DIR/Male_R2_ab.kmer21 -threads $THREADS

      meryl -B -C -m 21 -s $OUTPUT_DIR/${FEMALE_READ_1}_aa -o $OUTPUT_DIR/Female_R1_aa.kmer21 -threads $THREADS
      meryl -B -C -m 21 -s $OUTPUT_DIR/${FEMALE_READ_1}_ab -o $OUTPUT_DIR/Female_R1_ab.kmer21 -threads $THREADS

      meryl -B -C -m 21 -s $OUTPUT_DIR/${FEMALE_READ_2}_aa -o $OUTPUT_DIR/Female_R2_aa.kmer21 -threads $THREADS
      meryl -B -C -m 21 -s $OUTPUT_DIR/${FEMALE_READ_2}_ab -o $OUTPUT_DIR/Female_R2_ab.kmer21 -threads $THREADS

      echo "Removing split fastq files"
      rm -f $OUTPUT_DIR/*_aa $OUTPUT_DIR/*_ab

      echo "Adding kmer tables from different fastq chunks"
      meryl -M add -s $OUTPUT_DIR/Male_R1_aa.kmer21 -s $OUTPUT_DIR/Male_R1_ab.kmer21 -o $OUTPUT_DIR/Male_R1.kmer21
      meryl -M add -s $OUTPUT_DIR/Male_R2_aa.kmer21 -s $OUTPUT_DIR/Male_R2_ab.kmer21 -o $OUTPUT_DIR/Male_R2.kmer21
      meryl -M add -s $OUTPUT_DIR/Female_R1_aa.kmer21 -s $OUTPUT_DIR/Female_R1_ab.kmer21 -o $OUTPUT_DIR/Female_R1.kmer21
      meryl -M add -s $OUTPUT_DIR/Female_R2_aa.kmer21 -s $OUTPUT_DIR/Female_R2_ab.kmer21 -o $OUTPUT_DIR/Female_R2.kmer21

      echo "Adding kmer tables from paired reads"
      meryl -M add -s $OUTPUT_DIR/Male_R1.kmer21 -s $OUTPUT_DIR/Male_R2.kmer21 -o $OUTPUT_DIR/Male.kmer21
      meryl -M add -s $OUTPUT_DIR/Female_R1.kmer21 -s $OUTPUT_DIR/Female_R2.kmer21 -o $OUTPUT_DIR/Female.kmer21

      echo "Removing low-occurence kmers"
      meryl -M greaterthanorequal 2 -s $OUTPUT_DIR/Male.kmer21 -o $OUTPUT_DIR/Male.kmer21.cut2
      meryl -M greaterthanorequal 2 -s $OUTPUT_DIR/Female.kmer21 -o $OUTPUT_DIR/Female.kmer21.cut2
      meryl -M and -s $OUTPUT_DIR/Male.kmer21.cut2 -s $OUTPUT_DIR/Female.kmer21.cut2 -o $OUTPUT_DIR/Common.kmer21.cut2

      echo "kmer-masking the genome"
      kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms 21 -edb $OUTPUT_DIR/Male.kmer21.mask -mdb $OUTPUT_DIR/Male.kmer21.cut2 -l $ROUNDED_MAX_LENGTH -1 $OUTPUT_DIR/$FAKE_FASTQ -o $OUTPUT_DIR/Male_cut2
      kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms 21 -edb $OUTPUT_DIR/Female.kmer21.mask -mdb $OUTPUT_DIR/Female.kmer21.cut2 -l $ROUNDED_MAX_LENGTH -1 $OUTPUT_DIR/$FAKE_FASTQ -o $OUTPUT_DIR/Female_cut2
      kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms 21 -edb $OUTPUT_DIR/Common.kmer21.mask -mdb $OUTPUT_DIR/Common.kmer21.cut2 -l $ROUNDED_MAX_LENGTH -1 $OUTPUT_DIR/$FAKE_FASTQ -o $OUTPUT_DIR/Common_cut2

      echo "Writing out results"
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Male_cut2.match.1.fastq > $OUTPUT_DIR/Male_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Male_cut2.clean.1.fastq >> $OUTPUT_DIR/Male_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Female_cut2.match.1.fastq > $OUTPUT_DIR/Female_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Female_cut2.clean.1.fastq >> $OUTPUT_DIR/Female_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Common_cut2.match.1.fastq > $OUTPUT_DIR/Common_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Common_cut2.clean.1.fastq >> $OUTPUT_DIR/Common_cut2.tsv

    else
      echo "Using existing kmer databases"

      echo "kmer-masking the genome"
      kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms 21 -edb $OUTPUT_DIR/Male.kmer21.mask -mdb $DB_DIR/Male.kmer21.cut2 -l $ROUNDED_MAX_LENGTH -1 $OUTPUT_DIR/$FAKE_FASTQ -o $OUTPUT_DIR/Male_cut2
      kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms 21 -edb $OUTPUT_DIR/Female.kmer21.mask -mdb $DB_DIR/Female.kmer21.cut2 -l $ROUNDED_MAX_LENGTH -1 $OUTPUT_DIR/$FAKE_FASTQ -o $OUTPUT_DIR/Female_cut2
      kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms 21 -edb $OUTPUT_DIR/Common.kmer21.mask -mdb $DB_DIR/Common.kmer21.cut2 -l $ROUNDED_MAX_LENGTH -1 $OUTPUT_DIR/$FAKE_FASTQ -o $OUTPUT_DIR/Common_cut2

      echo "Writing out results"
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Male_cut2.match.1.fastq > $OUTPUT_DIR/Male_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Male_cut2.clean.1.fastq >> $OUTPUT_DIR/Male_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Female_cut2.match.1.fastq > $OUTPUT_DIR/Female_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Female_cut2.clean.1.fastq >> $OUTPUT_DIR/Female_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Common_cut2.match.1.fastq > $OUTPUT_DIR/Common_cut2.tsv
      kmer_match_genome_sizes.pl Genome_sizes.tsv $OUTPUT_DIR/Common_cut2.clean.1.fastq >> $OUTPUT_DIR/Common_cut2.tsv
    fi

elif  [ -z "$MALE_READ_1" ] || [ -z "$MALE_READ_2" ] || [ -z "$FEMALE_READ_1" ] || [ -z "$FEMALE_READ_2" ] || [ -z "$GENOME_INPUT" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required options."
    print_help
fi
