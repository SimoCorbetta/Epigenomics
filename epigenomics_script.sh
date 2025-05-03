#!/bin/bash

rep1=""
rep2=""
control=""
encode=""
black_list=""
qvalue=""
local_flag=""
filtered_flag=""
macs2_flag=""

while [ "$#" -gt 0 ]; do
    case $1 in
        -q|--qvalue)
            qvalue="$2"
            shift 2
            ;;
        -l|--local)
            local_flag="--local"
            shift
            ;;
        --filtered)
            filtered_flag="--filtered"
            shift
            ;;
        --macs2)
            macs2_flag="--macs2"
            shift
            ;;
        *)
            break
            ;;
    esac
done

if [ "$local_flag" ]; then
    # If --local option is present, only rep1, rep2, and control are required
    rep1=$1
    rep2=$2
    control=$3

    # Check if rep1, rep2, and control are provided
    if [ -z "$rep1" ] || [ -z "$rep2" ] || [ -z "$control" ]; then
        echo "Usage with local files: $0 [-l|--local] [--filtered] [--macs2] ENC_rep1 ENC_rep2 ENC_control"
        echo "It expects to find ENC_rep1.bam, ENC_rep2.bam, ENC_control.bam, black_list.bed and encode.bed in the pwd"
        exit 1
    fi
elif [ "$filtered_flag" ]; then
    echo "Using local filtered files (--filtered flag)"
    echo "It expects to find rep1_filter.bam, rep2_filter.bam, control_filter.bam, black_list.bed and encode.bed in the pwd"
    rep1="rep1"
    rep2="rep2"
    control="control"
elif [ "$macs2_flag" ]; then
    echo "Using local files outputted by MACS2 (--macs2 flag)"
    rep1="rep1"
    rep2="rep2"
    control="control"
else
    rep1=$1
    rep2=$2
    control=$3
    encode=$4
    black_list=$5

    # Check if all parameters are provided when not using --local
    if [ -z "$rep1" ] || [ -z "$rep2" ] || [ -z "$control" ] || [ -z "$encode" ] || [ -z "$black_list" ]; then
        echo "Download files usage: $0 [-q|--qvalue threshold] [-l|--local] [--filtered] [--macs2] ENC_rep1 ENC_rep2 ENC_control ENC_encode ENC_black_list"
        exit 1
    fi
fi


if [ "$qvalue" ]; then
    # If -q or --qvalue option is present, do something with the quality value
    echo -e "\e[32mq-value threshold specified: ${qvalue}\e[0m"
else
    qvalue=0.01
    echo -e "\e[32mUsing default q-value threshold: ${qvalue}\e[0m"
fi

if [ "$local_flag" ] || [ "$filtered_flag" ] || [ "$macs2_flag" ]; then
    # If --local option is present, use local files and rename them
    echo -e "\e[32mUsing local files...\e[0m"

else
    # Download files from ENCODE
    echo -e "\e[32mDownloading files...\e[0m"
    # wget rep1:
    wget "https://www.encodeproject.org/files/${rep1}/@@download/${rep1}.bam" -O ${rep1}.bam
    # wget rep2:
    wget "https://www.encodeproject.org/files/${rep2}/@@download/${rep2}.bam" -O ${rep2}.bam
    # wget control:
    wget "https://www.encodeproject.org/files/${control}/@@download/${control}.bam" -O ${control}.bam
    # wget for blacklist:
    wget "https://www.encodeproject.org/files/${black_list}/@@download/${black_list}.bed.gz" -O black_list.bed.gz
    # wget for encode peaks:
    wget "https://www.encodeproject.org/files/${encode}/@@download/${encode}.bed.gz" -O encode_not_sorted.bed.gz
    # wget for chromHMM:
    wget "http://159.149.160.56/HepG2_chromHMM_18states_hg38.bed" -O HepG2_chromHMM_18states_hg38.bed
fi

cp "${rep1}.bam" rep1.bam
cp "${rep2}.bam" rep2.bam
cp "${control}.bam" control.bam

gunzip *.gz

# Sort encode.bed
sort -k1,1 -k2,2n encode_not_sorted.bed> encode.bed


stat_file=stats_${rep1}_${rep2}_${control}.txt


# Check if either --filtered_flag or --macs2_flag is not set
if [ -z "${filtered_flag}" ] && [ -z "${macs2_flag}" ]; then

    # Do something when --filtered_flag is not set
    # Filtration of peaks

    # Filter peaks for replicates and control
    echo -e "\e[32mFiltering peaks for rep1\e[0m"
    samtools view -bq 1 rep1.bam > rep1_filter.bam
    echo -e "\e[32mFiltering peaks for rep2\e[0m"
    samtools view -bq 1 rep2.bam > rep2_filter.bam
    echo -e "\e[32mFiltering peaks for control\e[0m"
    samtools view -bq 1 control.bam > control_filter.bam

    # Contrast reads after filtering

    # Create a stats file
    echo -e "\e[34mGenerating statistics in ${stat_file}\e[0m"
    echo "statistics" > ${stat_file}

    # Statistics for the first replicate
    echo "Stats for rep1" >> ${stat_file}
    echo "number of reads before filtering"  >> ${stat_file}
    samtools view -c rep1.bam >> ${stat_file}
    echo "number of reads after filtering"  >> ${stat_file}
    samtools view -c rep1_filter.bam >> ${stat_file}

    # Statistics for the second replicate
    echo -e "Stats for rep2" >> ${stat_file}
    echo "number of reads before filtering"  >> ${stat_file}
    samtools view -c rep2.bam >> ${stat_file}
    echo "number of reads after filtering"  >> ${stat_file}
    samtools view -c rep2_filter.bam >> ${stat_file}

    # Statistics for the control
    echo -e "Stats for control" >> ${stat_file}
    echo "number of reads before filtering"  >> ${stat_file}
    samtools view -c control.bam >> ${stat_file}
    echo "number of reads after filtering"  >> ${stat_file}
    samtools view -c control_filter.bam >> ${stat_file}
fi

# Check if --macs2_flag is not set
if [ -z "${macs2_flag}" ]; then
    # Do something when --macs2_flag is not set
    # Peak calling

    # Call peaks using MACS2 for each replicate and merged
    echo -e "\e[33mCalling peaks for rep1\e[0m"
    macs2 callpeak -t rep1_filter.bam -c control_filter.bam -q ${qvalue} -g hs -n rep1 -B --SPMR > /dev/null 2> /dev/null
    echo -e "\e[33mCalling peaks for rep2\e[0m"
    macs2 callpeak -t rep2_filter.bam -c control_filter.bam -q ${qvalue} -g hs -n rep2 -B --SPMR > /dev/null 2> /dev/null
    echo -e "\e[33mCalling peaks for merged\e[0m"
    macs2 callpeak -t rep1_filter.bam rep2_filter.bam -c control_filter.bam -q ${qvalue} -g hs -n merged -B --SPMR > /dev/null 2> /dev/null

    # Regular expressions to filter only canonical chromosomes (1-22, X, Y)
    canonical_chromosomes='\b^chr([1-9]|1[0-9]|2[0-2]|X|Y)\b'

    # Filtering peaks and summits files
    echo -e "\e[32mFiltering peaks files\e[0m"
    grep -E "${canonical_chromosomes}" rep1_peaks.narrowPeak > rep1_peaks_chr_filt.narrowPeak
    grep -E "${canonical_chromosomes}" rep2_peaks.narrowPeak > rep2_peaks_chr_filt.narrowPeak
    grep -E "${canonical_chromosomes}" merged_peaks.narrowPeak > merged_peaks_chr_filt.narrowPeak

    echo -e "\e[32mFiltering summits files\e[0m"
    grep -E "${canonical_chromosomes}" rep1_summits.bed > rep1_summits_chr_filt.bed
    grep -E "${canonical_chromosomes}" rep2_summits.bed > rep2_summits_chr_filt.bed
    grep -E "${canonical_chromosomes}" merged_summits.bed > merged_summits_chr_filt.bed

     # Filtering treat pileup files
    echo -e "\e[32mFiltering treat pileup files\e[0m"
    grep -E "${canonical_chromosomes}" merged_treat_pileup.bdg > merged_treat_chr_filt_pileup.bdg
    grep -E "${canonical_chromosomes}" rep1_treat_pileup.bdg > rep1_treat_chr_filt_pileup.bdg
    grep -E "${canonical_chromosomes}" rep2_treat_pileup.bdg > rep2_treat_chr_filt_pileup.bdg

fi


# Statistics on peaks for the first replicate
echo -e "\e[34mGenerating statistics for rep1 peaks\e[0m"
echo "statistics filtering non-canonical chromosomes" >> ${stat_file}
echo "rep1" >> ${stat_file}
echo "number of reads before filtering" >> ${stat_file}
wc -l rep1_peaks.narrowPeak >> ${stat_file}
echo "number of reads after filtering"  >> ${stat_file}
wc -l rep1_peaks_chr_filt.narrowPeak >> ${stat_file}

# Statistics on peaks for the second replicate
echo -e "\e[34mGenerating statistics for rep2 peaks\e[0m"
echo "rep2" >> ${stat_file}
echo "number of reads before filtering" >> ${stat_file}
wc -l rep2_peaks.narrowPeak >> ${stat_file}
echo "number of reads after filtering"  >> ${stat_file}
wc -l rep2_peaks_chr_filt.narrowPeak >> ${stat_file}

# Statistics on merged peaks
echo -e "\e[34mGenerating statistics for merged peaks\e[0m"
echo "merged" >> ${stat_file}
echo "number of reads before filtering" >> ${stat_file}
wc -l merged_peaks.narrowPeak >> ${stat_file}
echo "number of reads after filtering"  >> ${stat_file}
wc -l merged_peaks_chr_filt.narrowPeak >> ${stat_file}

# Removing unfiltered files
echo -e "\e[31mRemoving unfiltered files\e[0m"
rm *_summits.bed
rm *_lambda.bdg
rm *_peaks.narrowPeak



# Use bedtools subtract to remove peaks in the black list
echo -e "\e[31mRemoving blacklisted peaks\e[0m"
bedtools subtract -a rep1_peaks_chr_filt.narrowPeak -b black_list.bed > rep1_peaks_white_listed.narrowPeak
bedtools subtract -a rep2_peaks_chr_filt.narrowPeak -b black_list.bed > rep2_peaks_white_listed.narrowPeak
bedtools subtract -a merged_peaks_chr_filt.narrowPeak -b black_list.bed > merged_peaks_white_listed.narrowPeak

# Obtain summits without blacklisted peaks
echo -e "\e[32mGenerating summit files without blacklisted peaks\e[0m"
awk '{print $1 "\t" $2+$10 "\t" $2+$10+1 "\t" $4 "\t" $7 "\t" $9'} rep1_peaks_white_listed.narrowPeak > rep1_summits_white_listed.bed
awk '{print $1 "\t" $2+$10 "\t" $2+$10+1 "\t" $4 "\t" $7 "\t" $9'} rep2_peaks_white_listed.narrowPeak > rep2_summits_white_listed.bed
awk '{print $1 "\t" $2+$10 "\t" $2+$10+1 "\t" $4 "\t" $7 "\t" $9'} merged_peaks_white_listed.narrowPeak > merged_summits_white_listed.bed

# Generate statistics on peaks after removing blacklisted regions
echo -e "\e[34mGenerating statistics after removing blacklisted peaks\e[0m"
echo "statistics filtering peaks in the black list" >> ${stat_file}

# Statistics for the first replicate
echo "rep1" >> ${stat_file}
echo "number of reads before filtering black list" >> ${stat_file}
wc -l rep1_peaks_chr_filt.narrowPeak >> ${stat_file}
echo "number of reads after filtering"  >> ${stat_file}
wc -l rep1_peaks_white_listed.narrowPeak >> ${stat_file}

# Statistics for the second replicate
echo "rep2" >> ${stat_file}
echo "number of reads before filtering black list" >> ${stat_file}
wc -l rep2_peaks_chr_filt.narrowPeak >> ${stat_file}
echo "number of reads after filtering"  >> ${stat_file}
wc -l rep2_peaks_white_listed.narrowPeak >> ${stat_file}

# Statistics for the merged peaks
echo "merged" >> ${stat_file}
echo "number of reads before filtering black list" >> ${stat_file}
wc -l merged_peaks_chr_filt.narrowPeak >> ${stat_file}
echo "number of reads after filtering"  >> ${stat_file}
wc -l merged_peaks_white_listed.narrowPeak >> ${stat_file}

# Remove blacklisted narrow peaks files
echo -e "\e[31mRemoving blacklisted narrow peaks files\e[0m"
rm *_chr_filt.narrowPeak



# Get the number of peaks in each replicate
rep1_peaks_number=$(wc -l < rep1_peaks_white_listed.narrowPeak)
rep2_peaks_number=$(wc -l < rep2_peaks_white_listed.narrowPeak)

# Determine which replicate has more peaks
if [ ${rep1_peaks_number} -ge ${rep2_peaks_number} ]; then
  rep_a=rep2_peaks_white_listed.narrowPeak
  rep_b=rep1_peaks_white_listed.narrowPeak
else
  rep_a=rep1_peaks_white_listed.narrowPeak
  rep_b=rep2_peaks_white_listed.narrowPeak
fi

# Perform bedtools intersect to find common peaks between replicates
echo -e "\e[32mFinding intersection between replicates\e[0m"
bedtools intersect -a ${rep_a} -b ${rep_b} -u > rep1rep2.narrowPeak

# Perform bedtools intersect and subtract to find peaks in rep1 not overlapping with encode
echo -e "\e[32mFinding intersection and subtracting peaks\e[0m"
bedtools intersect -a rep1_peaks_white_listed.narrowPeak -b encode.bed -u > rep1encode.narrowPeak
bedtools subtract -a rep1_peaks_white_listed.narrowPeak -b rep1encode.narrowPeak > rep1not_inencode.narrowPeak

# Perform bedtools intersect and subtract to find peaks in rep2 not overlapping with encode
bedtools intersect -a rep2_peaks_white_listed.narrowPeak -b encode.bed -u > rep2encode.narrowPeak
bedtools subtract -a rep2_peaks_white_listed.narrowPeak -b rep2encode.narrowPeak > rep2not_inencode.narrowPeak

# Perform bedtools intersect and subtract to find peaks in merged not overlapping with encode
bedtools intersect -a merged_peaks_white_listed.narrowPeak -b encode.bed -u > mergeencode.narrowPeak
bedtools subtract -a merged_peaks_white_listed.narrowPeak -b mergeencode.narrowPeak > mergenot_inencode.narrowPeak

# Perform bedtools intersect and subtract to find peaks in rep1rep2 not overlapping with encode
bedtools intersect -a rep1rep2.narrowPeak -b encode.bed -u > rep1rep2encode.narrowPeak
bedtools subtract -a rep1rep2.narrowPeak -b rep1rep2encode.narrowPeak > rep1rep2not_inencode.narrowPeak

# Statistics on the number of common peaks for each intersection
echo -e "\e[34mGenerating statistics on intersection peaks\e[0m"
echo "peaks found in intersect rep1 and rep2" >> ${stat_file}
wc -l rep1rep2.narrowPeak >> ${stat_file}
echo "peaks found in intersect rep1 and encode" >> ${stat_file}
wc -l rep1encode.narrowPeak >> ${stat_file}
echo "peaks found in intersect rep2 and encode" >> ${stat_file}
wc -l rep2encode.narrowPeak >> ${stat_file}
echo "peaks found in intersect merged and encode" >> ${stat_file}
wc -l mergeencode.narrowPeak >> ${stat_file}
echo "peaks found in intersect rep1_rep2 (already intersected) and encode" >> ${stat_file}
wc -l rep1rep2encode.narrowPeak >> ${stat_file}


# Find peaks in rep1 whose summit is within 100bp of a summit in rep2
echo -e "\e[32mFinding peaks in summit proximity between rep1 and rep2\e[0m"
bedtools closest -a rep1_summits_white_listed.bed -b rep2_summits_white_listed.bed -d -t first > rep1_rep2_closest_summits.bed
awk '{ if ($13 < 100 && $13 >= 0) { print } }' rep1_rep2_closest_summits.bed > rep1_rep2_closest_summits100.bed

# Find peaks in rep1 whose summit is within 100bp of encode summits
echo -e "\e[32mFinding peaks in summit proximity between rep1 and encode\e[0m"
awk '{ print $1 "\t" $2 + $10 "\t" $2 + $10 + 1}' encode.bed > encode_summits_not_sorted.bed
sort -k1,1 -k2,2n encode_summits_not_sorted.bed> encode_summits.bed
bedtools closest -a rep1_summits_white_listed.bed -b encode_summits.bed -d -t first > rep1_encode_closest_summits.bed
awk '{ if ($10 < 100 && $10 >= 0) { print } }' rep1_encode_closest_summits.bed > rep1_encode_closest_summits100.bed
bedtools subtract -a rep1_summits_white_listed.bed -b rep1_encode_closest_summits100.bed > rep1_not_encode_closest_summits100.bed

# Find peaks in rep2 whose summit is within 100bp of encode summits
echo -e "\e[32mFinding peaks in summit proximity between rep2 and encode\e[0m"
bedtools closest -a rep2_summits_white_listed.bed -b encode_summits.bed -d -t first > rep2_encode_closest_summits.bed
awk '{ if ($10 < 100 && $10 >= 0) { print } }' rep2_encode_closest_summits.bed > rep2_encode_closest_summits100.bed
bedtools subtract -a rep2_summits_white_listed.bed -b rep2_encode_closest_summits100.bed > rep2_not_encode_closest_summits100.bed

# Find peaks in merged whose summit is within 100bp of encode summits
echo -e "\e[32mFinding peaks in summit proximity between merged and encode\e[0m"
bedtools closest -a merged_summits_white_listed.bed -b encode_summits.bed -d -t first > merged_encode_closest_summits.bed
awk '{ if ($10 < 100 && $10 >= 0) { print } }' merged_encode_closest_summits.bed > merged_encode_closest_summits100.bed
bedtools subtract -a merged_summits_white_listed.bed -b merged_encode_closest_summits100.bed > merged_not_encode_closest_summits100.bed

# Find peaks in rep1rep2 whose summit is within 100bp of encode summits
echo -e "\e[32mFinding peaks in summit proximity between rep1rep2 and encode\e[0m"
bedtools closest -a rep1rep2.narrowPeak -b encode.bed -d -t first > rep1rep2_encode_closest_summits.bed
awk '{ if ($21 < 100 && $21 >= 0) { print } }' rep1rep2_encode_closest_summits.bed > rep1rep2_encode_closest_summits100.bed
bedtools subtract -a rep1rep2.narrowPeak -b rep1rep2_encode_closest_summits100.bed > rep1rep2_not_encode_closest_summits100.bed

# Statistics on the number of peaks in common for each summit proximity intersection
echo -e "\e[34mGenerating statistics on summit proximity peaks\e[0m"
echo "peaks found in summit proximity rep1 and rep2" >> ${stat_file}
wc -l rep1_rep2_closest_summits100.bed >> ${stat_file}
echo "peaks found in summit proximity rep1 and encode" >> ${stat_file}
wc -l rep1_encode_closest_summits100.bed >> ${stat_file}
echo "peaks found in summit proximity rep2 and encode" >> ${stat_file}
wc -l rep2_encode_closest_summits100.bed >> ${stat_file}
echo "peaks found in summit proximity merged and encode" >> ${stat_file}
wc -l merged_encode_closest_summits100.bed >> ${stat_file}
echo "peaks found in summit proximity rep1_rep2 and encode" >> ${stat_file}
wc -l rep1rep2_encode_closest_summits100.bed >> ${stat_file}



echo -e "\e[32mCutting pvalue and qvalue for peak overlapping\e[0m"
cut -f 7,9 rep1rep2.narrowPeak > peak_over_rep1rep2_p_q.tsv
cut -f 7,9 rep1encode.narrowPeak > peak_over_rep1encode_p_q.tsv
cut -f 7,9 rep2encode.narrowPeak > peak_over_rep2encode_p_q.tsv
cut -f 7,9 mergeencode.narrowPeak > peak_over_mergeencode_p_q.tsv
cut -f 7,9 rep1rep2encode.narrowPeak > peak_over_rep1rep2encode_p_q.tsv

echo -e "\e[32mCutting pvalue and qvalue for replicates not intersecting with encode (peak overlapping)\e[0m"
cut -f 7,9 rep1not_inencode.narrowPeak > rep1not_inencode.tsv
cut -f 7,9 rep2not_inencode.narrowPeak > rep2not_inencode.tsv
cut -f 7,9 mergenot_inencode.narrowPeak > mergenot_inencode.tsv
cut -f 7,9 rep1rep2not_inencode.narrowPeak > rep1rep2not_inencode.tsv

echo -e "\e[32mCutting pvalue and qvalue for summit proximity\e[0m"
cut -f 5,6 rep1_rep2_closest_summits100.bed > sum_prox_rep1rep2_p_q.tsv
cut -f 5,6 rep1_encode_closest_summits100.bed > sum_prox_rep1encode_p_q.tsv
cut -f 5,6 rep2_encode_closest_summits100.bed > sum_prox_rep2encode_p_q.tsv
cut -f 5,6 merged_encode_closest_summits100.bed > sum_prox_mergeencode_p_q.tsv
cut -f 7,9 rep1rep2_encode_closest_summits100.bed > sum_prox_rep1rep2encode_p_q.tsv

echo -e "\e[32mCutting pvalue and qvalue for replicates not intersecting with encode (summit proximity)\e[0m"
cut -f 5,6 rep1_not_encode_closest_summits100.bed > rep1not_inencode_summits100.tsv
cut -f 5,6 rep2_not_encode_closest_summits100.bed > rep2not_inencode_summits100.tsv
cut -f 5,6 merged_not_encode_closest_summits100.bed > mergenot_inencode_summits100.tsv
cut -f 7,9 rep1rep2not_inencode.narrowPeak > rep1rep2not_inencode.tsv

echo -e "\e[32mCreating boxplots for comparisons between files\e[0m"
Rscript "boxplotting_new.R" peak_over_rep1encode_p_q.tsv rep1not_inencode.tsv "rep1 vs encode" > /dev/null 2> /dev/null
Rscript "boxplotting_new.R" peak_over_rep2encode_p_q.tsv rep2not_inencode.tsv "rep2 vs encode" > /dev/null 2> /dev/null
Rscript "boxplotting_new.R" peak_over_mergeencode_p_q.tsv mergenot_inencode.tsv "merged vs encode" > /dev/null 2> /dev/null
Rscript "boxplotting_new.R" peak_over_rep1rep2encode_p_q.tsv rep1rep2not_inencode.tsv "rep1rep2 vs encode" > /dev/null 2> /dev/null

Rscript "boxplotting_new.R" sum_prox_rep1encode_p_q.tsv rep1not_inencode_summits100.tsv "rep1 vs encode (summit proximity)" > /dev/null 2> /dev/null
Rscript "boxplotting_new.R" sum_prox_rep2encode_p_q.tsv rep2not_inencode_summits100.tsv "rep2 vs encode (summit proximity)" > /dev/null 2> /dev/null
Rscript "boxplotting_new.R" sum_prox_mergeencode_p_q.tsv mergenot_inencode_summits100.tsv "merged vs encode (summit proximity)" > /dev/null 2> /dev/null
Rscript "boxplotting_new.R" sum_prox_rep1rep2encode_p_q.tsv rep1rep2not_inencode.tsv "rep1rep2 vs encode (summit proximity)" > /dev/null 2> /dev/null



echo -e "\e[32mSorting BED files for Jaccard index calculation\e[0m"
sort -k1,1 -k2,2n encode.bed > encode_peaks_sorted.bed
sort -k1,1 -k2,2n merged_peaks_white_listed.narrowPeak > merged_peaks_white_listed_sorted.narrowPeak
sort -k1,1 -k2,2n rep1rep2.narrowPeak > rep1rep2_sorted.narrowPeak

echo -e "\e[32mCalculating Jaccard Index for merged peaks and encode peaks\e[0m"
result_merged=$(bedtools jaccard -a merged_peaks_white_listed_sorted.narrowPeak -b encode_peaks_sorted.bed)
# Extract the Jaccard index value
jaccard_index_merged=$(echo "${result_merged}" | awk 'NR==2 {print $3}')
echo "Jaccard Index for merged peaks and encode peaks: ${jaccard_index_merged}"

echo -e "\e[32mCalculating Jaccard Index for rep1rep2 peaks and encode peaks\e[0m"
result_intersect=$(bedtools jaccard -a rep1rep2_sorted.narrowPeak -b encode_peaks_sorted.bed)
# Extract the Jaccard index value
jaccard_index_intersect=$(echo "${result_intersect}" | awk 'NR==2 {print $3}')
echo "Jaccard Index for rep1rep2 peaks and encode peaks: ${jaccard_index_intersect}"

echo -e "\e[32mChoosing the peaks file with the higher Jaccard Index\e[0m"
if [ ${jaccard_index_merged} -ge ${jaccard_index_intersect} ]; then
  rep_a=merged_peaks_white_listed_sorted.narrowPeak
else
  rep_a=rep1rep2_sorted.narrowPeak
fi
echo "Selected peaks file: ${rep_a}"

echo -e "\e[32mCalculating intersection with chromatin state annotations\e[0m"
bedtools intersect -a ${rep_a} -b HepG2_chromHMM_18states_hg38.bed -wa -wb > chromatin_states.bed

echo -e "\e[32mGenerating barplot using Python script\e[0m"
# Run Python script to generate barplot
python3 barplot_pie.py chromatin_states.bed

echo -e "\e[32mGenerating file for GREAT tool\e[0m"
# Create a narrowPeak file for GREAT tool (using only the first 4 columns)
cut -f 1-4 ${rep_a} > def_peaks_forGreat.narrowPeak
