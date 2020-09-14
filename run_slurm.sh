#!/usr/bin/bash
#SBATCH --output=logs/G34-%j.out
#SBATCH --time=00:10:00
#SBATCH --account=rrg-jabado-ab
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

mkdir -p logs/

module load r/3.5.1

# usage:
# sbatch --export=ALL,rmd='01-my_analysis.Rmd' run_slurm.sh
R --no-save -e "rmarkdown::render('${rmd}', 'html_document')"
