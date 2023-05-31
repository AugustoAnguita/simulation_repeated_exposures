#! /bin/bash
#SBATCH --job-name=s1_exp3 
#SBATCH --partition=long
#SBATCH --account=generic
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=augusto.anguita@isglobal.org 
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=24
#SBATCH --mem=96gb
#SBATCH --output=/PROJECTES/ATHLETE/ATHLETE_sim_long/DLNM/ANALYSIS/simulation_mars_2022/dataY1andX/exp3/logs/dataY1exp3_extramodels.log # Standard output and error log

module purge > /dev/null 2>&1
module load lang/R/4.0.3-foss-2020a

Rscript data1exp3_new_extramodels.R
