#!/usr/bin/bash
#PBS -N anderson
#PBS -o logs/
#PBS -e logs/
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=50G
#PBS -l vmem=50G
#PBS -l epilogue=/mnt/KLEINMAN_BACKUP/home/selin.jessa/epilogue.sh
 
cd $PBS_O_WORKDIR
mkdir -p logs/

export TMPDIR="/mnt/KLEINMAN_BACKUP/home/selin.jessa/tmp"
export R_LIBS_USER="/mnt/KLEINMAN_JBOD1/KLEINMAN_BACKUP/home/selin.jessa/renv_global/renv/library/R-3.5/x86_64-redhat-linux-gnu"

# usage: qsub -v rmd=02-explore_joining.Rmd run_hydra.sh
R --no-save -e "rmarkdown::render('${rmd}', 'html_document')"
