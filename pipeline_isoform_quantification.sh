#Step 1: Counting isoform specific kmer from the fastq.gz trimmed data

for i in $(cat /home/data/qs/HGC20231021002-0004/trimmed_data/sample_list.txt)
do
  echo "Processing sample: $i"
  python /home/data/qs/HGC20231021002-0004/trimmed_data/kmer_count_opti_v2.py  --fastq_path /home/data/qs/HGC20231021002-0004/trimmed_data/${i}_trim_1P.fastq.gz --num_threads 20 --chunk_size 100000 --csv_input_dir /home/data/qs/data/reference_isoform/mus_ref/50mer_Mus_musculus_GRCc39_isoform_mini_shared_filter --csv_output_dir /home/data/qs/HGC20231021002-0004/trimmed_data/output_csv_${i}
done

# Step 2: Merging the isoform specific reference csv and corresponding counting csv file together for downstream analysis

for i in $(cat /home/data/qs/HGC20231021002-0004/trimmed_data/sample_list.txt)
do
  echo "Merging files for sample: $i"
  python /home/data/qs/HGC20231021002-0004/trimmed_data/merge_normalizing_count.py --kmer_reference_directory /home/data/qs/data/reference_isoform/mus_ref/50mer_Mus_musculus_GRCc39_isoform_mini_shared_filter --kmer_counts_directory /home/data/qs/HGC20231021002-0004/trimmed_data/output_csv_${i} --fastq /home/data/qs/HGC20231021002-0004/trimmed_data/${i}_trim_1P.fastq.gz --output_directory /home/data/qs/HGC20231021002-0004/trimmed_data/${i}_merged_data --read_length 150 --k 50
done

# Step 3: GC based GaussF for isoform quantificaiton
for i in $(cat /home/data/qs/HGC20231021002-0004/trimmed_data/sample_list.txt)
do
  echo "isoform quantification: $i"
  python /home/data/qs/HGC20231021002-0004/trimmed_data/GaussF_esti.py --threshold 20 --input /home/data/qs/HGC20231021002-0004/trimmed_data/${i}_merged_data/ --output /home/data/qs/HGC20231021002-0004/trimmed_data/results_file_GC_${i}.csv
done

