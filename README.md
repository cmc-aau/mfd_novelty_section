# mfd_novelty_section
This repo is part of the [Microflora Danica](https://github.com/cmc-aau/mfd_wiki/wiki) (MFD) project and it reproduces the results from "Section VI: MAGs from Danish habitats double the known species fraction in metagenomes" of the [MFD manuscripts](https://www.biorxiv.org/content/10.1101/2024.06.27.600767v1).

## Input data
Please download the input file from [Zenodo](). 

The SingleM metapackakge with GTDB R214 supplemented with other genomes can be downloaded from this [Zenodo](https://zenodo.org/records/12741285).

The SingleM metapackakge with GTDB R214 supplemented with other genomes, together with our short-read (SR) MFD MAGs can be downloaded from this [Zenodo]().

## Codes and scripts
All SingleM condensed profiles from our big batch samples are merged together using the script: <code>/scripts/scripts_R/00-merge_condensed_singlem_outputs.r</code>. 

R scripts used for generating SingleM-related plots (eg., Figure 5a, Figure 5b, Supplementary Figure 6, and Supplementary Figure 7) in the paper are in the folder <code>/scripts/scripts_R</code>.

## SingleM example commands: 
Taking "sample1" and metapackage "uhgg_smag_spire_oceans_gems.smpkg" as an example:
### 1. SingleM microbial profiling: SingleM pipe metagenome
```
singlem pipe -1 /PATH_TO_FORWARD_READ/sample1_trimmed.fastq.gz -2 /PATH_TO_REVERSE_READ/sample1_trimmed.fastq.gz --otu-table sample1.otu_table.csv --threads 10 --metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg --archive-otu-table sample1.otu_table.csv.archive --taxonomic-profile sample1.condensed_otu_defaultcutoffs
```
### 2. SingleM microbial fraction approximation
```
singlem microbial_fraction --input-profile sample1.condensed_otu_defaultcutoffs --output-tsv sample1_readfrac_stdout.tsv --output-per-taxon-read-fractions sample1_readfrac_pertaxon -1 /PATH_TO_FORWARD_READ/sample1_trimmed.fastq.gz -2 /PATH_TO_REVERSE_READ/sample1_trimmed.fastq.gz --metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
```
### 3. SingleM for MAG recovery of microbial community approximation
Taking "MAG1" as an example:
#### 3a. SingleM pipe MAGs
```
singlem pipe --genome-fasta-files /PATH_TO_FILE/MAG1.fa --otu-table MAG1_genome.otu_table.csv --threads 2 --archive-otu-table MAG1_genome.otu_table_archive.csv --taxonomic-profile MAG1_condensed_otu_defaultcutoffs --metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
```
#### 3b. SingleM summarise to combine otu tables of all MAGs
```
singlem summarise --input-otu-tables-list /PATH_TO_FILE/ABSOLUTE_FILEPATHS_TO_MAG_OTUTABLES_LIST.txt --output-otu-table genomes.otu_table_combined.csv --metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
```
#### 3b. SingleM summarise to combine otu tables of all metagenomes
```
singlem summarise --input-otu-tables-list /PATH_TO_FILE/ABSOLUTE_FILEPATHS_TO_METAGENOME_OTUTABLES_LIST.txt --output-otu-table metagenome_otu_table_combinedotu.txt --metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg
```
#### 3c. SingleM appraise
```
singlem appraise --metagenome-otu-tables /PATH_TO_FILE/metagenome_otu_table_combinedotu.txt --genome-otu-tables /PATH_TO_FILE/genomes.otu_table_combined.csv --metapackage /PATH_TO_FILE/uhgg_smag_spire_oceans_gems.smpkg --output-binned-otu-table appraise_binned_otu_table.csv --output-unbinned-otu-table appraise_unbinned_otu_table.csv --output-unaccounted-for-otu-table appraise_unaccounted_for.csv --output-found-in &> singlem_appraise.log &
```
