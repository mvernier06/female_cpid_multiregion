#! /bin/bash
#SBATCH --time=0-24:00:00
#SBATCH -p public
#SBATCH --output=run_R_job%J.out
#SBATCH --error=run_R_job%J.err
#SBATCH --mail-type=END  
#SBATCH --mail-user=marinevernier@unistra.fr

# SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G

# Configuration de l'environnement
source $(conda info --base)/etc/profile.d/conda.sh
conda activate r_4.4

checktime=$(date);
echo "start time $checktime";

# Run the R script
Rscript /home2020/home/inci/mvernier/cpid_multireg_female/female_cpid_multiregion/4__MEGENA/5_modules_enrichment.R 
checktime=$(date);
echo "end time $checktime"; 

exit 0