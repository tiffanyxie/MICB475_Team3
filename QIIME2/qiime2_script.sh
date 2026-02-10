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

# After visualization, decided to trim at 390
# Determine ASVs with DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
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

# Train classifier for V1-V3
#use truncal length of 390 used for trimming
#Use universal 27F (5′- AGA GTT TGA TCM TGG CTC AG–3′) and
# 519R (5′- GWA TTA CCG CGG CKG CTG–3′) primers
qiime feature-classifier extract-reads \
  --i-sequences /datasets/classifiers/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer AGAGTTTGATCMTGGCTCAG \
  --p-r-primer GWATTACCGCGGCKGCTG \
  --p-trunc-len 390 \
  --o-reads ref-seqs-trimmed.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-trimmed.qza \
  --i-reference-taxonomy /datasets/classifiers/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier.qza