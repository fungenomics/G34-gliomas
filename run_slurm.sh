#!/usr/bin/bash
#SBATCH --job-name="G34-gliomas"
#SBATCH --time=01:00:00
#SBATCH --output=logs/myproject-%j.out
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G

mkdir -p logs/

module load r/3.5.1

# usage:
# sbatch --export=ALL,rmd='01-my_analysis.Rmd' run_slurm.sh
R --no-save -e "rmarkdown::render('${rmd}', 'html_document')"
# Rscript singlecell_tumor/analysis/02-get_sample_table_V2.R
