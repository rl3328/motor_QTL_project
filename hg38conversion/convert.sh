#!/bin/bash

# # download tranform file
# cd /mnt/lustre/home/rl3328/rl3328/motor_qtl/GWAS

# # download LiftOver tools
# wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
# chmod +x liftOver

# # download hg19 → hg38 chain file
# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz


##############  main step to do transformation from hg19 to hg38
# original hg19 annotation
INPUT_FILE="./file_path.txt"

# outputs
OUTPUT_PATH="/mnt/lustre/home/rl3328/rl3328/motor_qtl/GWAS"
UNMAPPED_PATH="${OUTPUT_PATH}/unmapped"

mkdir -p "${OUTPUT_PATH}"
mkdir -p "${UNMAPPED_PATH}"

# read input txt, for each path in this txt, do liftOver transformation
while IFS= read -r file_path; do
    if [[ -f "$file_path" ]]; then
        # extract basename（remove .bed）
        base_name=$(basename "$file_path" .bed)

        #  liftOver
        ./liftOver \
            "$file_path" \
            /mnt/lustre/home/rl3328/rl3328/motor_qtl/GWAS/hg19ToHg38.over.chain.gz \
            "${OUTPUT_PATH}/${base_name}.to_hg38.bed" \
            "${UNMAPPED_PATH}/${base_name}.to_hg38.unmapped.bed"

        echo "Converted: $file_path → ${OUTPUT_PATH}/${base_name}.to_hg38.bed"
    else
        echo "File not found: $file_path"
    fi
done < "$INPUT_FILE"

echo "Conversion completed!"