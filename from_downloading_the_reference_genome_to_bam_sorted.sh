# O link para download do hg38 está disponível em https://support.illumina.com/sequencing/sequencing_software/igenome.html
%%time
!wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
!tar -xf Homo_sapiens_UCSC_hg38.tar.gz

# codigo index
bowtie2-build referencia.fasta nome_do_indice_bowtie2
bwa index referencia.fasta
bwa-mem2 index referencia.fasta


# download from dbgap
while read sra_id; do prefetch --ngc "your_project.ngc" "$sra_id"; done < "arquivo_unico.txt"

# fastq dump
fastq-dump --gzip --split-files SAMPLE3.sra

# alinhamento com o bwa-mem 2 e .sam to .bam
bwa-mem2 mem /path/to/workdir/genome.fa SAMPLE1_1.fastq.gz SAMPLE1_2.fastq.gz | samtools view -b -o SAMPLE1_bwa2.bam

# alinhamento com o bwa-mem e .sam to .bam
bwa mem /path/to/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa SAMPLE4_1.fastq.gz SAMPLE4_2.fastq.gz | samtools view -bS -1 - > SAMPLE4_mapped.bam

# alinhamento com o Bowtie2 e .sam to .bam
bowtie2 -x /path/to/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.fa -1 SAMPLE4_1.fastq.gz -2 SAMPLE4_2.fastq.gz | samtools view -b -o SAMPLE4_bowtie.bam

# sort e index com o samtools
samtools sort SAMPLE3_bwa2.bam -o SAMPLE3_bwa2_sorted.bam

samtools index SAMPLE1_bwa2_sorted.bam
