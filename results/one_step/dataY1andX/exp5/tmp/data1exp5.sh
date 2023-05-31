#! /bin/bash
#SBATCH --job-name=ATHLETE_longsim_s1_exp5
#SBATCH --partition=long
#SBATCH --account=generic
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=augusto.anguita@isglobal.org 
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=30
#SBATCH --mem=140gb
#SBATCH --output=/PROJECTES/SHARED/AugustoAnguita/ATHLETE_sim_long/DLNM/ANALYSIS/simulation_mars_2022/dataY1andX/exp5/logs/dataY1exp5.log # Standard output and error log

module purge > /dev/null 2>&1
module load lang/R/4.0.3-foss-2020a

Rscript data1exp5_new.R
