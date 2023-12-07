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
    echo "  -kmer_length             Set the kmer size (default: 21))"
    echo "  -h                       Print this help message"
    exit 0
}

SKIP_COUNTING=false
DB_DIR=""
SKIP_SPLITTING=false
THREADS=1
KMER_LENGTH=21


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
        -kmer_length)
            shift
            KMER_LENGTH=$1
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

if [ "$SKIP_COUNTING" = true ] && [ -f "$DB_DIR/Male.kmer${KMER_LENGTH}.cut2.mcdat" ] && [ -f "$DB_DIR/Female.kmer${KMER_LENGTH}.cut2.mcdat" ] && [ -f "$DB_DIR/Common.kmer${KMER_LENGTH}.cut2.mcdat" ] && [ -f "$GENOME_INPUT" ] && [ ! -z "$OUTPUT_DIR" ]; then
    echo "No read files needed, reading from existing kmer databases"
    CHECK=true

else
  echo "mcdat files do not exist"
    CHECK=false
fi

if [  "$SKIP_COUNTING" = false ] && [ -f "$MALE_READ_2" ] && [ -f "$FEMALE_READ_1" ] && [ -f "$FEMALE_READ_2" ] && [ -f "$GENOME_INPUT" ] && [ ! -z "$OUTPUT_DIR" ]; then
	echo "Data files ok, will generate new kmer databases"
  CHECK=true

else
  echo "Error: One or more input files do not exist."
  print_help
  CHECK=false
fi

if [ "$CHECK" = true ]; then

      if [ ! -d "$OUTPUT_DIR" ]; then
          mkdir -p "$OUTPUT_DIR" || { echo "Error: Unable to create output directory $OUTPUT_DIR"; exit 1; }
      fi


      ln -s $GENOME_INPUT $OUTPUT_DIR/Genome_file.fasta
      GENOME_INPUT=$OUTPUT_DIR/Genome_file.fasta


      if [ "$SKIP_SPLITTING" = false ]; then

            echo "Splitting genome to Ns"
            split_genome_to_Ns.pl $GENOME_INPUT
            GENOME="${GENOME_INPUT}.split"
            GENOME_INDEX="${GENOME}.fai"
            FAKE_FASTQ=$(echo $GENOME | sed 's/\.split$/.fastq/')


      else
        echo "Genome will not be split"
        GENOME=$GENOME_INPUT
        GENOME_INDEX="${GENOME}.fai"
        FAKE_FASTQ=$(echo $GENOME | sed 's/\.fasta$/.fastq/')
      fi


      echo "Generating fake fastq genome file"
      fasta2fake_fastq.pl $GENOME

      echo "Indexing the genome"
      samtools faidx $GENOME


      cut -f1,2 $GENOME_INDEX > $OUTPUT_DIR/Genome_sizes.tsv

       echo "Calculating max sequence length"
          MAX_LENGTH=$(awk -F'\t' -v col="2" '{ if ($col > max) max = $col } END { print max }' "$OUTPUT_DIR/Genome_sizes.tsv")
          ROUNDED_MAX_LENGTH=$(echo "scale=0; ((($MAX_LENGTH + 99) / 100) * 100)" | bc)


      if [ "$SKIP_COUNTING" = false ]; then
        ln -s $MALE_READ_1 $OUTPUT_DIR/male_read_1
        ln -s $MALE_READ_2 $OUTPUT_DIR/male_read_2
        ln -s $FEMALE_READ_1 $OUTPUT_DIR/female_read_1
        ln -s $FEMALE_READ_2 $OUTPUT_DIR/female_read_2

        MALE_READ_1=$OUTPUT_DIR/male_read_1
        MALE_READ_2=$OUTPUT_DIR/male_read_2
        FEMALE_READ_1=$OUTPUT_DIR/female_read_1
        FEMALE_READ_1=$OUTPUT_DIR/female_read_1


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


          split -l $LINES_PER_FILE_MALE_1 $MALE_READ_1 ${MALE_READ_1}_
          split -l $LINES_PER_FILE_MALE_2 $MALE_READ_2 ${MALE_READ_2}_
          split -l $LINES_PER_FILE_FEMALE_1 $FEMALE_READ_1 ${FEMALE_READ_1}_
          split -l $LINES_PER_FILE_FEMALE_2 $FEMALE_READ_2 ${FEMALE_READ_2}_

          echo "Generating kmer table from reads"
          meryl -B -C -m $KMER_LENGTH -s ${MALE_READ_1}_aa -o $OUTPUT_DIR/Male_R1_aa.kmer${KMER_LENGTH} -threads $THREADS
          meryl -B -C -m $KMER_LENGTH -s ${MALE_READ_1}_ab -o $OUTPUT_DIR/Male_R1_ab.kmer${KMER_LENGTH} -threads $THREADS

          meryl -B -C -m $KMER_LENGTH -s ${MALE_READ_2}_aa -o $OUTPUT_DIR/Male_R2_aa.kmer${KMER_LENGTH} -threads $THREADS
          meryl -B -C -m $KMER_LENGTH -s ${MALE_READ_2}_ab -o $OUTPUT_DIR/Male_R2_ab.kmer${KMER_LENGTH} -threads $THREADS

          meryl -B -C -m $KMER_LENGTH -s ${FEMALE_READ_1}_aa -o $OUTPUT_DIR/Female_R1_aa.kmer${KMER_LENGTH} -threads $THREADS
          meryl -B -C -m $KMER_LENGTH -s ${FEMALE_READ_1}_ab -o $OUTPUT_DIR/Female_R1_ab.kmer${KMER_LENGTH} -threads $THREADS

          meryl -B -C -m $KMER_LENGTH -s ${FEMALE_READ_2}_aa -o $OUTPUT_DIR/Female_R2_aa.kmer${KMER_LENGTH} -threads $THREADS
          meryl -B -C -m $KMER_LENGTH -s ${FEMALE_READ_2}_ab -o $OUTPUT_DIR/Female_R2_ab.kmer${KMER_LENGTH} -threads $THREADS

          echo "Removing split fastq files"
          rm -f $OUTPUT_DIR/*_aa $OUTPUT_DIR/*_ab

          echo "Adding kmer tables from different fastq chunks"
          meryl -M add -s $OUTPUT_DIR/Male_R1_aa.kmer${KMER_LENGTH} -s $OUTPUT_DIR/Male_R1_ab.kmer${KMER_LENGTH} -o $OUTPUT_DIR/Male_R1.kmer${KMER_LENGTH}
          meryl -M add -s $OUTPUT_DIR/Male_R2_aa.kmer${KMER_LENGTH} -s $OUTPUT_DIR/Male_R2_ab.kmer${KMER_LENGTH} -o $OUTPUT_DIR/Male_R2.kmer${KMER_LENGTH}
          meryl -M add -s $OUTPUT_DIR/Female_R1_aa.kmer${KMER_LENGTH} -s $OUTPUT_DIR/Female_R1_ab.kmer${KMER_LENGTH} -o $OUTPUT_DIR/Female_R1.kmer${KMER_LENGTH}
          meryl -M add -s $OUTPUT_DIR/Female_R2_aa.kmer${KMER_LENGTH} -s $OUTPUT_DIR/Female_R2_ab.kmer${KMER_LENGTH} -o $OUTPUT_DIR/Female_R2.kmer${KMER_LENGTH}

          echo "Adding kmer tables from paired reads"
          meryl -M add -s $OUTPUT_DIR/Male_R1.kmer${KMER_LENGTH} -s $OUTPUT_DIR/Male_R2.kmer${KMER_LENGTH} -o $OUTPUT_DIR/Male.kmer${KMER_LENGTH}
          meryl -M add -s $OUTPUT_DIR/Female_R1.kmer${KMER_LENGTH} -s $OUTPUT_DIR/Female_R2.kmer${KMER_LENGTH} -o $OUTPUT_DIR/Female.kmer${KMER_LENGTH}

          echo "Removing low-occurence kmers"
          meryl -M greaterthanorequal 2 -s $OUTPUT_DIR/Male.kmer${KMER_LENGTH} -o $OUTPUT_DIR/Male.kmer${KMER_LENGTH}.cut2
          meryl -M greaterthanorequal 2 -s $OUTPUT_DIR/Female.kmer${KMER_LENGTH} -o $OUTPUT_DIR/Female.kmer${KMER_LENGTH}.cut2
          meryl -M and -s $OUTPUT_DIR/Male.kmer${KMER_LENGTH}.cut2 -s $OUTPUT_DIR/Female.kmer${KMER_LENGTH}.cut2 -o $OUTPUT_DIR/Common.kmer${KMER_LENGTH}.cut2

          echo "kmer-masking the genome"
          kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms $KMER_LENGTH -edb $OUTPUT_DIR/Male.kmer${KMER_LENGTH}.mask -mdb $OUTPUT_DIR/Male.kmer${KMER_LENGTH}.cut2 -l $ROUNDED_MAX_LENGTH -1 $FAKE_FASTQ -o $OUTPUT_DIR/Male_cut2
          kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms $KMER_LENGTH -edb $OUTPUT_DIR/Female.kmer${KMER_LENGTH}.mask -mdb $OUTPUT_DIR/Female.kmer${KMER_LENGTH}.cut2 -l $ROUNDED_MAX_LENGTH -1 $FAKE_FASTQ -o $OUTPUT_DIR/Female_cut2
          kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms $KMER_LENGTH -edb $OUTPUT_DIR/Common.kmer${KMER_LENGTH}.mask -mdb $OUTPUT_DIR/Common.kmer${KMER_LENGTH}.cut2 -l $ROUNDED_MAX_LENGTH -1 $FAKE_FASTQ -o $OUTPUT_DIR/Common_cut2

          echo "Writing out results"
          kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Male_cut2.match.1.fastq > $OUTPUT_DIR/Male_cut2.tsv
          kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Male_cut2.clean.1.fastq >> $OUTPUT_DIR/Male_cut2.tsv
          kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Female_cut2.match.1.fastq > $OUTPUT_DIR/Female_cut2.tsv
          kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Female_cut2.clean.1.fastq >> $OUTPUT_DIR/Female_cut2.tsv
          kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Common_cut2.match.1.fastq > $OUTPUT_DIR/Common_cut2.tsv
          kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Common_cut2.clean.1.fastq >> $OUTPUT_DIR/Common_cut2.tsv

      else
        echo "Using existing kmer databases"

        echo "kmer-masking the genome"
        kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms $KMER_LENGTH -edb $OUTPUT_DIR/Male.kmer${KMER_LENGTH}.mask -mdb $DB_DIR/Male.kmer${KMER_LENGTH}.cut2 -l $ROUNDED_MAX_LENGTH -1 $FAKE_FASTQ -o $OUTPUT_DIR/Male_cut2
        kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms $KMER_LENGTH -edb $OUTPUT_DIR/Female.kmer${KMER_LENGTH}.mask -mdb $DB_DIR/Female.kmer${KMER_LENGTH}.cut2 -l $ROUNDED_MAX_LENGTH -1 $FAKE_FASTQ -o $OUTPUT_DIR/Female_cut2
        kmer-mask -nomasking -t $THREADS -v -match 0.2 -clean 0.2 -e 10 -m 3 -ms $KMER_LENGTH -edb $OUTPUT_DIR/Common.kmer${KMER_LENGTH}.mask -mdb $DB_DIR/Common.kmer${KMER_LENGTH}.cut2 -l $ROUNDED_MAX_LENGTH -1 $FAKE_FASTQ -o $OUTPUT_DIR/Common_cut2

        echo "Writing out results"
        kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Male_cut2.match.1.fastq > $OUTPUT_DIR/Male_cut2.tsv
        kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Male_cut2.clean.1.fastq >> $OUTPUT_DIR/Male_cut2.tsv
        kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Female_cut2.match.1.fastq > $OUTPUT_DIR/Female_cut2.tsv
        kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Female_cut2.clean.1.fastq >> $OUTPUT_DIR/Female_cut2.tsv
        kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Common_cut2.match.1.fastq > $OUTPUT_DIR/Common_cut2.tsv
        kmer_match_genome_sizes.pl $OUTPUT_DIR/Genome_sizes.tsv $OUTPUT_DIR/Common_cut2.clean.1.fastq >> $OUTPUT_DIR/Common_cut2.tsv
      fi

elif  [ "$CHECK"=false ]; then
    echo "Error: Missing required options."
    print_help
fi
