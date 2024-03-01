#!/bin/bash

echo "Processing bacterial samples"
BATCH=ume_hb_bac_samples
# Set up primers
FWD=GTGYCAGCMGCCGCGGTAA
REV=GGACTACNVGGGTWTCTAAT
RC_FWD=TTACCGCGGCKGCTGRCAC
RC_REV=ATTAGAWACCCBNGTAGTCC

# Create file with sample names to loop over
ls *{forward,R1_001}.fastq.gz | cut -f1 -d "." | awk -F'_R1_001' '{print $1}' > ${BATCH}.txt

## Choose the end of where the primers will be trimmed
## Documentation available at https://cutadapt.readthedocs.io/en/stable/guide.html#regular-3-adapters
export FIVE_PRIME=g
export THREE_PRIME=a
export FIVE_PRIME_PAIR=G
export THREE_PRIME_PAIR=A

## For multi-threading or parallel computing set a variable for number of threads
NCORES=${2:-1}

# Trim the reads
for sample in $(cat $BATCH.txt)
do

  echo "On sample: $sample"
  if test -f ${sample}.forward.fastq.gz; then
    FILE_PAT_F=.forward
    FILE_PAT_R=.reverse
  else 
    FILE_PAT_F=_R1_001
    FILE_PAT_R=_R2_001
  fi  
  F=${sample}${FILE_PAT_F}.fastq.gz
  R=${sample}${FILE_PAT_R}.fastq.gz
  cutadapt -${THREE_PRIME} ${FWD}...${RC_REV} \
    -${THREE_PRIME_PAIR} ${REV}...${RC_FWD} \
    -m 150 -M 302 -j ${NCORES}\
    -o ${sample}_R1_trimmed.fq.gz -p ${sample}_R2_trimmed.fq.gz \
    ${F} ${R} \
    >> cutadapt_${BATCH}_trimming_stats.txt 2>&1
done

##*** Clean up the folder a bit
mkdir trimmed_fastq_${BATCH}
mv *trimmed.fq.gz trimmed_fastq_${BATCH}/

##*** Check out how much you kept after trimming
##* This first line creates column titles for text file
##* The rest of the script will extract the sample name, the number of read pairs 
##* written, the percent of reads that passed filtering, and the percent of reads kept after filtering
echo $'Parameters: minimum of 150bp and maximum of 300 bp\nSample\t\tpairs_processed\ttoo_short\ttoo_long\tpairs_written\tpassing_filters\tfiltered' > summary_${BATCH}_filtered.txt
paste ${BATCH}.txt <(grep "improperly paired" cutadapt_${BATCH}_trimming_stats.txt | cut -f2 -d ".") <(grep "read pairs" cutadapt_${BATCH}_trimming_stats.txt | cut -f2 -d ":") <(grep "too short" cutadapt_${BATCH}_trimming_stats.txt | cut -f2 -d ":" | cut -f1 -d "(") <(grep "too long" cutadapt_${BATCH}_trimming_stats.txt | cut -f2 -d ":" | cut -f1 -d "(") <(grep "passing" cutadapt_${BATCH}_trimming_stats.txt | cut -f2 -d ":" | cut -f1 -d "(") <(grep "passing" cutadapt_${BATCH}_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_${BATCH}_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") > summary_${BATCH}_filtered.append.txt
cat summary_${BATCH}_filtered.txt summary_${BATCH}_filtered.append.txt > summary_${BATCH}_filter.txt
rm summary_${BATCH}_filtered.append.txt summary_${BATCH}_filtered.txt

# ## To find which ones were imporperly paired (need to find a better way to code this)
paste <(grep "improperly paired" cutadapt_${BATCH}_trimming_stats.txt | cut -f2 -d ".") > ${BATCH}_improperly_paired.txt

##*** Now that you're done with the original fastqs put them in a file
mkdir demultiplexed_og_${BATCH}_fastqs
mv *.fastq.gz demultiplexed_og_${BATCH}_fastqs/
