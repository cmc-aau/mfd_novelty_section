# mfd_novelty_section
This repo is part of the [Microflora Danica](https://github.com/cmc-aau/mfd_wiki/wiki) (MFD) project and it reproduces the results from "Section VI: MAGs from Danish habitats double the known species fraction in metagenomes" of the [MFD manuscripts](https://www.biorxiv.org/content/10.1101/2024.06.27.600767v1).

## Input data
Please download the input file from [Zenodo](). 

The SingleM metapackakge with GTDB R214 supplemented with other genomes can be downloaded from this [Zenodo](https://zenodo.org/records/12741285).

The SingleM metapackakge with GTDB R214 supplemented with other genomes, together with our short-read (SR) MFD MAGs can be downloaded from this [Zenodo]().

## Codes and scripts
All SingleM condensed profiles from our big batch samples are merged together using the script: <code>/scripts/00-merge_condensed_singlem_outputs.r</code>. 

R scripts used for generating SingleM-related plots (eg., Figure 5a, Figure 5b, Supplementary Figure 6, and Supplementary Figure 7) in the paper are in the folder <code>/scripts</code>.

## SingleM example commands: 
Taking "sample1" and metapackage "uhgg_smag_spire_oceans_gems.smpkg" as an example:
### 1. SingleM microbial profiling: SingleM pipe metagenome
```
singlem pipe
-1 /PATH_TO_FORWARD_READ/sample1_trimmed.fastq.gz
-2 /PATH_TO_REVERSE_READ/sample1_trimmed.fastq.gz
--threads 10
--metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
--otu-table sample1.otu_table.csv
--archive-otu-table sample1.otu_table.csv.archive
--taxonomic-profile sample1.condensed_otu_defaultcutoffs
```
### 2. SingleM microbial fraction approximation
```
singlem microbial_fraction
--input-profile sample1.condensed_otu_defaultcutoffs
--output-tsv sample1_readfrac_stdout.tsv
--output-per-taxon-read-fractions sample1_readfrac_pertaxon
-1 /PATH_TO_FORWARD_READ/sample1_trimmed.fastq.gz
-2 /PATH_TO_REVERSE_READ/sample1_trimmed.fastq.gz
--metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
```
### 3. SingleM for MAG recovery of microbial community approximation
Taking "MAG1" as an example:
#### 3a. SingleM pipe MAGs
```
singlem pipe
--genome-fasta-files /PATH_TO_FILE/MAG1.fa
--otu-table MAG1_genome.otu_table.csv
--threads 2
--archive-otu-table MAG1_genome.otu_table_archive.csv
--taxonomic-profile MAG1_condensed_otu_defaultcutoffs
--metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
```
#### 3b. SingleM summarise to combine otu tables of all MAGs
```
singlem summarise
--input-otu-tables-list /PATH_TO_FILE/ABSOLUTE_FILEPATHS_TO_MAG_OTUTABLES_LIST.txt
--output-otu-table genomes.otu_table_combined.csv
--metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
```
#### 3b. SingleM summarise to combine otu tables of all metagenomes
```
singlem summarise
--input-otu-tables-list /PATH_TO_FILE/ABSOLUTE_FILEPATHS_TO_METAGENOME_OTUTABLES_LIST.txt
--output-otu-table metagenome_otu_table_combinedotu.txt
--metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
```
#### 3c. SingleM appraise
```
singlem appraise
--metagenome-otu-tables /PATH_TO_FILE/metagenome_otu_table_combinedotu.txt
--genome-otu-tables /PATH_TO_FILE/genomes.otu_table_combined.csv
--metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
--output-binned-otu-table appraise_binned_otu_table.csv
--output-unbinned-otu-table appraise_unbinned_otu_table.csv
--output-unaccounted-for-otu-table appraise_unaccounted_for.csv
--output-found-in &> singlem_appraise.log &
```
## Additional note:
SingleM results on NCBI datasets come from [Ben J. Woodcroft, et al. (2024) 's work](https://doi.org/10.1101/2024.01.30.578060):

SingleM and Sandpiper: Robust microbial taxonomic profiles from metagenomic data., Ben J. Woodcroft, Samuel T. N. Aroney, Rossen Zhao, Mitchell Cunningham, Joshua A. M. Mitchell, Linda Blackall, Gene W. Tyson., bioRxiv 2024.01.30.578060; doi: https://doi.org/10.1101/2024.01.30.578060

The <code>microbial_fraction</code> module of SingleM:

Eisenhofer, Raphael, Antton Alberdi, and Ben J. Woodcroft. Large-scale estimation of bacterial and archaeal DNA prevalence in metagenomes reveals biome-specific patterns. bioRxiv (2024): 2024-05. https://doi.org/10.1101/2024.05.16.594470

## Citation of MFD manuscript:
Microflora Danica: the atlas of Danish environmental microbiomes., CM Singleton, TBN Jensen, F Delogu, EA Sørensen, VR Jørgensen, SM Karst, Y Yang, KS Knudsen, M Sereika, F Petriglieri, S Knutsson, SM Dall, RH Kirkegaard, JM Kristensen, BJ Woodcroft, DR Speth, STN Aroney, The Microflora Danica Consortium, M Wagner, MKD Dueholm, PH Nielsen, M Albertsen., bioRxiv 2024.06.27.600767; doi: https://doi.org/10.1101/2024.06.27.600767
