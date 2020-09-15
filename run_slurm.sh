#!/usr/bin/bash
#SBATCH --time=01:00:00
#SBATCH --output=logs/myproject-%j.out
#SBATCH --account=rrg-jabado-ab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

mkdir -p logs/

module load r/3.5.1

# usage:
# sbatch --export=ALL,rmd='01-my_analysis.Rmd' run_slurm.sh
R --no-save -e "rmarkdown::render('${rmd}', 'html_document')"
