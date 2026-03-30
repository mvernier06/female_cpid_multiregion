#! /bin/bash
#SBATCH --time=0-24:00:00
# SBATCH -p public
#SBATCH -p grant -A g2025a419c
#SBATCH --output=run_R_job%J.out
#SBATCH --error=run_R_job%J.err
#SBATCH --mail-type=END  
#SBATCH --mail-user=marinevernier@unistra.fr
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --array=1-6

source $(conda info --base)/etc/profile.d/conda.sh
conda activate r_megena

checktime=$(date);
echo "start time $checktime";

regions=("ACC" "Ins" "Hb" "Nac")

pairs=(
"ACC Ins"
"ACC Hb"
"ACC Nac"
"Ins Hb"
"Ins Nac"
"Hb Nac"
)

pair=${pairs[$SLURM_ARRAY_TASK_ID-1]}

reg1=$(echo $pair | cut -d' ' -f1)
reg2=$(echo $pair | cut -d' ' -f2)

echo "Running $reg1 vs $reg2"

Rscript /home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/4__MEGENA/MEGENA_with_RIN_correction/module_preservation_between_regions/9_module_preservation_between_regions.R $reg1 $reg2
checktime=$(date);
echo "end time $checktime"; 

exit 0