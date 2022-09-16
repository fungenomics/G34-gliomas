#!/usr/bin/bash

# Selin Jessa
# April 26, 2018
# From https://github.com/czi-hca-comp-tools/easy-data/blob/master/datasets/ucsc_human_cortex.md

# Samples: We performed single-cell mRNA sequencing (scRNA-seq) in primary cortical and medial
# ganglionic eminence (MGE) samples across stages of peak neurogenesis.
# To capture cells along the span of neuronal differentiation from progenitor cells to postmitotic
# neurons, we included samples of microdissected germinal zone and cortical plate.
# In addition, to compare cells of the same type across distinct cortical areas, we captured cells from
# the prefrontal cortex (PFC) and primary visual cortex (V1), including 13 paired specimens.
# Library: 7137 cells were isolated with Fluidigm C1s and sequenced on Illumina HiSeq. 4261 cells
# were retained after filtered for mitochondrial and ribosomal expression.
# Clustering: To resolve finer distinctions among cells and to relate gene expression variation to
# cell type, developmental stage, and anatomical source, we performed weighted gene coexpression network
# analysis (WGCNA) (fig. S3 and tables S6 and S7) (8).
# We reasoned that although individual genes may not be reliably detected by using scRNA-seq because
# of high dropout rate, gene coexpression networks would be more robustly detected and could more
# accurately reflect the molecular variation related to distinct biological processes.

# Gene matrix values are TPMs from HISAT2, not log2-transformed.

# Note March 25th, 2019: the below links no longer work, see the following page for download
# https://cells.ucsc.edu/
wget http://cells.ucsc.edu/aparna/geneMatrix.tsv.gz
wget http://cells.ucsc.edu/aparna/meta.tsv
gunzip geneMatrix.tsv.gz
