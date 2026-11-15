#!/bin/bash

# Parte 1: Processamento com bcftools
#echo "Iniciando processamento com bcftools..."
#time (bcftools mpileup -f genome.fa SAMPLE1_bwa2_sorted.bam | bcftools call -mv | bgzip > SAMPLE1_bcftools.vcf.gz)
#time (bcftools mpileup -f genome.fa SAMPLE2_bwa2_sorted.bam | bcftools call -mv | bgzip > SAMPLE2_bcftools.vcf.gz)
#time (bcftools mpileup -f genome.fa SAMPLE3_bwa2_sorted.bam | bcftools call -mv | bgzip > SAMPLE3_bcftools.vcf.gz)

# Parte 2: Processamento com freebayes
echo "Iniciando processamento com freebayes..."
time (freebayes -f genome.fa SAMPLE1_bwa2_sorted.bam | bgzip > SAMPLE1_bwa2_sorted_freebayes.vcf.gz)
time (freebayes -f genome.fa SAMPLE2_bwa2_sorted.bam | bgzip > SAMPLE2_bwa2_sorted_freebayes.vcf.gz)
time (freebayes -f genome.fa SAMPLE3_bwa2_sorted.bam | bgzip > SAMPLE3_bwa2_sorted_freebayes.vcf.gz)

echo "Processamento concluído!"
