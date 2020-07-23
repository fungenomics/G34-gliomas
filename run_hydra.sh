#!/usr/bin/bash
#PBS -N G34
#PBS -o logs/
#PBS -e logs/
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=40G
#PBS -l vmem=40G
#PBS -l epilogue=/mnt/KLEINMAN_BACKUP/home/selin.jessa/epilogue.sh
 
cd $PBS_O_WORKDIR
mkdir -p logs/

# usage:
# $ qsub -V -N 13 -v script=13-integrate_across_tech_v3.R run_hydra_r_3.6.sh
# Rscript ${script}

# usage: qsub -v rmd=02-explore_joining.Rmd run_hydra.sh
R --no-save -e "rmarkdown::render('${rmd}', 'html_document')"
