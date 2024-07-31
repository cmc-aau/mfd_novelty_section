# mfd_novelty_section
This repo is part of the [Microflora Danica](https://github.com/cmc-aau/mfd_wiki/wiki) (MFD) project and it reproduces the results from "Section VI: MAGs from Danish habitats double the known species fraction in metagenomes" of the [MFD manuscripts](https://www.biorxiv.org/content/10.1101/2024.06.27.600767v1).

## Input data
Please download the input data from [Zenodo](). 
## Codes and scripts
All SingleM condensed profiles from our big batch samples are merged together using the script: <code>/scripts/scripts_R/00-merge_condensed_singlem_outputs.r</code>. 

R scripts used for generating SingleM-related plots (eg., Figure 5a, Figure 5b, Supplementary Figure 6, and Supplementary Figure 7) in the paper are in the folder <code>/scripts/scripts_R</code>.

Bash codes for running SingleM are provided in the folder <code>/scripts/scripts_bash</code>. Please note that SingleM scripts provided are meant to be run with a SLURM system and the SingleM program is downloaded through mamba.


