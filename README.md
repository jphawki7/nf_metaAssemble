# nf_metaAssemble
Nextflow pipeline for nanopore metagenomic assemblies and read classification\
Author: **Justin Paul Hawkins**

This pipeline will default perform adapter and quality trimming using **Porechop** and **Fastplong**, read classification using **Kraken2**, metagenomic De Novo assemblies using **Flye, Canu, and MetaMDBG**, and identify assembled contigs using **DIAMOND**.

Options also exsist for read de-hosting with **Hostile**, read filtering by taxid with **Krakentools**, and **BlastN**. See help for full list of commands.

# Prerequisites
- This pipeline requires Conda for use. Environment file can be found in the envs folder.
- To use Kraken2, Diamond, or Blast, you will require a local database which you can indicate in the .nf file, or with the appropriate flag
- If running kraken2 ensure your computers has enough RAM, pipeline is set to only run 1 kraken2 analysis at a time
- Ensure computer specs are reflected within the .conf file

# Quick Start
To run default analysis, activate conda environment, change directory to path with .nf, and run:
```
nextflow run nf_metaAssemble.nf --analysis <abs path to folder wtih .fastq reads> \
--krakendb <abs to kraken2 db> \
--diamonddb <abs to diamond db> \
-c config/base.config
```
# Flags
**--taxid <taxid>**        (Default 1, all reads, change to desired taxid to filter, requires kraken2)\
**--single**               (Set to only get requested taxid reads, skips unclassified and root reads\
**--skip_kraken**         (Skip kraken2 classification and read filtering)\
**--skip_diamond**         (Skip diamond search, useful if only blastn desired)\
**--skip_chop**           (Skip porechop read trimming, useful if reads already trimmed)\
**--hostile <abs path>**   (Abs path to host sequence for Hostile host removal - enables hostile filtering)\
**--outdir <abs path>**     (Abs path to output directory, default is analysis_output in read folder)\
**--blastn**             (Enable blastn search, requires --blastdb indicated, default diamond only)\
**--krakendb <abs path>**  (Abs path to kraken2 database)\
**--blastdb <abs path>**   (Abs path to blast database)\
**--diamonddb <abs path>** (Abs path to diamond database)\
