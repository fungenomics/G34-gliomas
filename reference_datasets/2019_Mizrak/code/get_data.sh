#!/usr/bin/bash

mkdir -p data
cd data

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109447/suppl/GSE109447_13055_Cell_Barcodes.xlsx
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109447/suppl/GSE109447_13055_cells_id.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109447/suppl/GSE109447_13055_cells_matrix.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109447/suppl/GSE109447_29319_Cell_Barcodes.xlsx
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109447/suppl/GSE109447_29319_cells.matrix.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109447/suppl/GSE109447_29319_cells_id_repinfo.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109447/suppl/GSE109447_Rep1_29319cells_Basic.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109447/suppl/GSE109447_Rep2_13055cells_Basic.txt.gz
