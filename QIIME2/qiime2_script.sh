#!/bin/bash

# Script for QIIME processing
# Last edited: Feb. 7, 2026 by Tiffany Xie

# Create directory for analysis outputs and navigate to it
mkdir /data/soil_project
cd /data/soil_project

#Import and demultiplex data using manifest file
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/soil/soil_manifest.tsv \
  --output-path ./demux_seqs.qza

#Create visualization of demultiplexed samples
qiime demux summarize \
  --i-data demux_seqs.qza \
  --o-visualization demux.qzv



