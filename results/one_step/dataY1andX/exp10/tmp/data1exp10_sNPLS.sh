#! /bin/bash
#SBATCH --job-name=ATHLETE_sNPLS_s1_exp10 
#SBATCH --partition=long
#SBATCH --account=generic
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=augusto.anguita@isglobal.org 
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=24
#SBATCH --mem=96gb
#SBATCH --output=/PROJECTES/ATHLETE/analyses/22_longi_meths_benchmark_AA/Simulation/COMPLETED_EXPERIMENT/MAIN_ANALYSIS/dataY1andX/exp10/logs/data1exp10_sNPLS.log # Standard output and error log

module purge > /dev/null 2>&1
module load lang/R/4.0.3-foss-2020a

Rscript data1exp10_sNPLS.R
