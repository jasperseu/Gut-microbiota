cat sample1_1.fastq.gz sample1_2.fastq.gz > sample1_cat.fastq.gz
..
cat samplen_1.fastq.gz samplen_2.fastq.gz > samplen_cat.fastq.gz


python "metaphlan2.py" "~/data/simple1_cat.fastq.gz" --input_type fastq --mpa_pkl "metaphlan_databases/mpa_v20_m200.pkl" --bowtie2db "metaphlan_databases/mpa_v20_m200" > "~/data/simple1_cat.taxa.txt"
..
python "metaphlan2.py" "~/data/simplen_cat.fastq.gz" --input_type fastq --mpa_pkl "metaphlan_databases/mpa_v20_m200.pkl" --bowtie2db "metaphlan_databases/mpa_v20_m200" > "~/data/simplen_cat.taxa.txt"

python merge_metaphlan_tables.py *taxa.txt > taxa_all.txt