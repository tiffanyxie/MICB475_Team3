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

# Use the trained classifier to assign taxonomy to your reads (rep-seqs.qza)
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

  #visualize taxonomy file
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

#Copy metadatafile into project folder
cp /datasets/project_2/soil/soil_metadata.txt .

#Taxa barplot
qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file soil_metadata.txt \
    --o-visualization taxa-bar-plots.qzv

#Filter out mitochondria and chloroplast ASVs
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file soil_metadata.txt

# Select for organic layer (o-horizon) and no herbicide treated (0) samples
qiime feature-table filter-samples \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file /datasets/project_2/soil/soil_metadata.txt \
  --p-where "[Horizon]='O horizon' AND [Herbicide Use]='0'" \
  --o-filtered-table o-layer-no-herbicide-filtered-table.qza

# Further filter for BC only table
  qiime feature-table filter-samples \
  --i-table o-layer-no-herbicide-filtered-table.qza \
  --m-metadata-file /datasets/project_2/soil/soil_metadata.txt \
  --p-where "[Region]='British Columbia'" \
  --o-filtered-table bc-only-o-layer-no-herbicide-filtered-table.qza

qiime feature-table summarize \
  --i-table bc-only-o-layer-no-herbicide-filtered-table.qza \
  --o-visualization bc-only-o-layer-no-herbicide-filtered-table.qzv \
  --m-sample-metadata-file /datasets/project_2/soil/soil_metadata.txt

# Generate tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

# only run until here (ivana) -----

# Alpha rarefaction curve (8000 is close to 8747)
qiime diversity alpha-rarefaction \
--i-table table.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 8000 \
--m-metadata-file /datasets/project_2/soil/soil_metadata.txt \
--o-visualization alpha-rarefaction.qzv

qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table.qza \
--p-sampling-depth 5296 \
--m-metadata-file /datasets/project_2/soil/soil_metadata.txt \
--output-dir core-metrics-results
