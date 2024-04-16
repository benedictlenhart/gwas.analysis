#!/bin/bash
#
###goal, creat a slurm script that gathers data from permutation files
#SBATCH -J residuals.perm # A single job name for the array
#SBATCH -c 10
#SBATCH -N 1 # on one node
#SBATCH -n 1 # one task
#SBATCH -t 3:30:00 #<= this may depend on your resources
#SBATCH --mem 20G #<= this may depend on your resources
#SBATCH -e /scratch/aob2x/score_error/gwas.gmmat.%A_%a.err # Standard error
#SBATCH -o /scratch/aob2x/score_output/gwas.gmmat.%A_%a.out # Standard output
#SBATCH -p standard
#SBATCH -A berglandlab_standard
### sbatch  --array=1 /standard/vol186/bergland-lab/Adam/gwas/gwas.analysis/gwas.nogrm.slurm.sh
### cat /scratch/aob2x/score_error/gwas.gmmat.
###9760 tasks tasks
### sacct -j 60223545


module load gcc/11.4
module load openmpi/4.1.4
module load R/4.3.1

#define variables
jobid=${SLURM_ARRAY_TASK_ID}
wd="/standard/vol186/bergland-lab/Adam/gwas/gwas.analysis/" #working directory
ncores=${SLURM_CPUS_PER_TASK}

Rscript ${wd}gwas.nogrm.perms.R  \
        ${jobid} \
        ${ncores}
