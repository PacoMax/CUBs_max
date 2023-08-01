# CUBs

Scripts for using R ENC'. A metric to quantify codon usage bias selection in prokaryote genomes. If you want to run them, please follow the instructions specified when running each script.

## CUBs_max.sh
Script for calculation of tAI and Nc' (metrics for CUBs (Codon usage bias selection)) for each protein-coding gene given a genome, and tables required by RENC_super.r.

## RENC_super.r
Script with functions for calculating R ENC', Codon usage bias selection (tAI and NC') on Gene Ontology categories, and also for doing GSEA for CUBs in GO categories. It uses outputs from CUBs_max.sh

## Scripts used by CUBs_max.sh

### get_table_protein.pl

### get_trnas_from_trnascan.pl

### prepare_tmp.pl


## Supplementary information for Master's Integrative Biology thesis dissertation

FigS1.svg

Figura S1.svg

Figura_S2.tif

Figura_S3.svg

Tabla S2.csv

Tabla_S1.csv


González Serrano, F. M. (2020). Evolución de la eficiencia traduccional debido a la selección del uso de codones en procariotas (Master's thesis, Tesis (MC)--Centro de Investigación y de Estudios Avanzados del IPN Unidad Irapuato. Departamento de Ingeniería Genética). https://repositorio.cinvestav.mx/bitstream/handle/cinvestav/1696/SSIT0016243.pdf?sequence=1

González-Serrano, F., Abreu-Goodger, C. & Delaye, L. Translation Comes First: Ancient and Convergent Selection of Codon Usage Bias Across Prokaryotic Genomes. J Mol Evol 90, 438–451 (2022). https://doi.org/10.1007/s00239-022-10074-0
