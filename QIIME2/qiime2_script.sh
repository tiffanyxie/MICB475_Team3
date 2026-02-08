#!/bin/bash

# Script for QIIME processing
# Last edited: Feb. 7, 2026 by Tiffany Xie & Ivana Djaja

# Create directory for analysis outputs and navigate to it
mkdir /data/soil_project
cd /data/soil_project

# Import and demultiplex data using manifest file
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/soil/soil_manifest.tsv \
  --output-path ./demux_seqs.qza

# Create visualization of demultiplexed samples
qiime demux summarize \
  --i-data demux_seqs.qza \
  --o-visualization demux.qzv

# ----- Ivana's add-on from here -----
# After visualization, decided to trim at 390
# Determine ASVs with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 390 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza

# Visualize ASVs stats
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file /datasets/project_2/soil/soil_metadata.txt

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv
