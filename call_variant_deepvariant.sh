#!/bin/bash
BIN_VERSION="1.6.1"
HOST_DIR="/path/to/workdir"  # Confirmed path

docker run \
  -v "${HOST_DIR}:/data" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="/data/genome.fa" \
  --reads="/data/SAMPLE1_bwa2_sorted.bam" \
  --output_vcf="/data/SAMPLE1_deepvariant.vcf.gz" \
  --num_shards=$(nproc)
