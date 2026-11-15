bcftools merge --force-samples -o var_deepvariants.vcf.gz -O z SAMPLE1_deepvariant.vcf.gz SAMPLE2_deepvariant.vcf.gz SAMPLE3_deepvariant.vcf.gz
bcftools merge --force-samples -o var_freebayes.vcf.gz -O z SAMPLE1_bwa2_sorted_freebayes.vcf.gz  SAMPLE2_bwa2_sorted_freebayes.vcf.gz  SAMPLE3_bwa2_sorted_freebayes.vcf.gz
bcftools merge -o var_bcftools.vcf.gz -O z SAMPLE1_bcftools.vcf.gz SAMPLE2_bcftools.vcf.gz SAMPLE3_bcftools.vcf.gz
bcftools merge -o var_gatk.vcf.gz -O z SAMPLE1_gatk.vcf.gz SAMPLE2_gatk.vcf.gz SAMPLE3_gatk.vcf.gz
