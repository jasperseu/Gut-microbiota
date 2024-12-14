cat *.fa > ref.fa
bowtie2-build ref.fa ref
for gz in `ls *.fastq.gz`;do name=`basename $gz .fastq.gz`;echo "bowtie2 -p 12 -x ../microbiome/ncbi_dataset/fna/ref --no-mixed --very-sensitive --n-ceil 0,0.01  ${name}.fastq.gz | samtools sort -O bam -@ 24 -o - >${name}ref.bam">>run_bowtie.command;done

for bam in `ls *bam`;do name=`basename $bam ref.bam`;echo "bcftools mpileup -C 50 -Ou -m 3 -F 0.0002 -f ./microbiome_T1/fna/ref.fa ${name}ref.bam > ${name}ref.bcf" >> run_bcftools1.command;done

for bam in `ls *bam`;do name=`basename $bam ref.bam`;echo "bcftools call -c --variants-only -Ob --ploidy 1 ${name}ref.bcf > ${name}.bcf" >> run_bcftools2.command;done_

for bam in `ls *bam`;do name=`basename $bam ref.bam`;echo "bcftools view -Ov  ${name}.bcf | vcfutils.pl varFilter -d 100 >  ${name}.vcf" >> run_bcftools3.command;done_

for bam in `ls *bam`;do name=`basename $bam ref.bam`;echo "bgzip ${name}.vcf" >> bgzip.command;done

for gz in `ls *.fastq.gz`;do name=`basename $gz .fastq.gz`;echo "java -Xmx10G -jar ~/jiayu/biosoft/snpEff/snpEff.jar eff -c ~/jiayu/biosoft/snpEff/snpEff.config AT_10 /media/desk16/sxq106/jiayu/RCC/rawdata\ T0/${name}.vcf.gz > ${name}.snp.eff.vcf -csvStats ${name}.csv -stats ${name}.html" >> run_snpeff.command;done