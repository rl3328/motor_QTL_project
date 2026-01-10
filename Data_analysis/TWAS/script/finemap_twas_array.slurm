#!/bin/bash
#SBATCH --job-name=twas_motor
#SBATCH --output=/mnt/lustre/home/rl3328/rl3328/TWAS/logs/mnm_array_%A_%a.out
#SBATCH --error=/mnt/lustre/home/rl3328/rl3328/TWAS/logs/mnm_array_%A_%a.err
#SBATCH --time=150:00:00
#SBATCH --mem=240GB
#SBATCH --array=0-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint="cpu2mem8a|cpu2mem8b|cpu2mem16a|cpu2mem16b|cpu4mem32a|cpu4mem32b|cpu8mem64a|cpu8mem64b|cpu32mem256a|cpu32mem256b"


export PATH="/mnt/lustre/lab/gwang/home/rl3328/.pixi/bin:$PATH"
source ~/.bashrc

TMPDIR=/mnt/lustre/lab/gwang/home/rl3328/tmp/slurm_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

home_dir=/mnt/lustre/lab/gwang/home/rl3328
cd ${home_dir}

trap "rm -rf $TMPDIR" EXIT

# Input file with all region info
REGION_FILE=/mnt/lustre/home/rl3328/rl3328/motor_qtl/twas_test/regions_missing_rds.txt
#/mnt/lustre/home/rl3328/rl3328/susie_ash/export/hg38_1362_blocks.bed
#/mnt/lustre/home/rl3328/rl3328/motor_qtl/regions_missing_rds.txt
#/mnt/lustre/home/rl3328/rl3328/motor_qtl/missing_regions.tsv
#
# cd ~/glyco_output/RSS_QC_nonimputed_result/twas && cut -f4 /mnt/lustre/home/rl3328/rl3328/susie_ash/export/hg38_1362_blocks.bed | while read b; do     clean_b=$(echo "$b" | tr -d '\r');     ls "ROSMAP_Bennett_Klein_glycoQTL.${clean_b}.twas.tsv.gz" &>/dev/null || echo "$clean_b"; done > ~/rl3328/gpQTL_project/TWAS/missing_regions.txt


# REGION_FILE=/mnt/lustre/home/rl3328/rl3328/gpQTL_project/TWAS/script/missing_regions.tsv

LINES_PER_JOB=1  # 302/31 ≈ 44, round up
START_LINE=$((SLURM_ARRAY_TASK_ID * LINES_PER_JOB + 1))
END_LINE=$((START_LINE + LINES_PER_JOB - 1))

for LINE_NUM in $(seq $START_LINE $END_LINE); do
    # Skip header line (line 1), so add 1 to LINE_NUM
    # REGION_NAME=$(awk -v line=$((LINE_NUM + 1)) 'NR==line {print $4}' $REGION_FILE) # depends on region file
    REGION_NAME=$(awk -v line=$((LINE_NUM + 1)) 'NR==line {print $1}' $REGION_FILE)
    [ -z "$REGION_NAME" ] && break  # Stop if past end of file

    echo "[$(date)] Processing region: $REGION_NAME"
    # Run your command
    sos run /mnt/lustre/lab/gwang/home/rl3328/xqtl-protocol/code/pecotmr_integration/twas_ctwas.ipynb twas \
        --name ROSMAP_eQTL_wp_unadjusted \
        --gwas_meta_data /mnt/lustre/home/rl3328/rl3328/motor_qtl/gwas_meta_unadjusted_wp_2020_raw.tsv \
        --ld_meta_data /mnt/lustre/data/gaowang/resource/LD_EUR/ld_meta_file.tsv \
        --xqtl_meta_data /mnt/lustre/home/rl3328/rl3328/motor_qtl/ROSMAP_eQTL.meta_twas.tsv \
        --xqtl_type_table /mnt/lustre/lab/gwang/home/rl3328/gpQTL_project/TWAS/data_type_table.txt \
        --rsq_pval_cutoff 0.05 --rsq_cutoff 0.01 --rsq_option adj_rsq \
        --cwd /mnt/lustre/lab/ctcn/hklein/motor_qtl_project/unadjusted_wp_2020/susie_twas/noQC_nonimputed_result/ -s build \
        --regions /mnt/lustre/lab/gwang/home/rl3328/resource/EUR_LD_blocks.bed \
        --region-name $REGION_NAME

done