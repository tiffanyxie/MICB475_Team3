# 1. QIIME2 Data Processing

## Aim:
* Import and demultiplex data from soil dataset
* Denoise data by trimming and removing low quality reads
* Cluster unique reads into amplicon sequence variants (ASVs)
* Filter out mitochondrial and chloroplast DNA

## Code
[QIIME2 Script](../QIIME2/QIIME2_script.sh)

## Results

### Post-demultiplexing
[demux.qzv](../QIIME2/output/demux.qzv)

**Demux Sequence Statistics**:

![Sequence counts summary](../QIIME2/output/demux_seq_counts.png)

**Quality Score per Base**:

![Base quality score](../QIIME2/output/demux_base_qual.png)

**Demux Sequence Lengths**:

![Sequence lengths](../QIIME2/output/demux_seq_length.png)

Decided to trim to 390 bp

### ASVs



<img width="870" height="590" alt="image" src="https://github.com/user-attachments/assets/744cd9fc-a4db-405b-87fc-26fec3fe3494" />

### Taxonomy Analysis
* Trained classifier using universal 27F and 519R primers covering V1 to V3
* Performed taxonomic analysis using this classifier
![Taxa bar plot](../QIIME2/output/level-4-bars.svg)

### Mitochondria and chloroplast filtered table\
* Lost 1 sample and ~ 7000 ASVs
<img width="699" height="555" alt="image" src="https://github.com/user-attachments/assets/30485b64-96a9-46a3-a114-fbb8b75d5353" />


## Next steps
* Decide what other metadata filtering needs to be done
* Perform rarefaction and generate diversity metrics

### Organic layer + no herbicide + BC only filtering
<img width="709" height="566" alt="image" src="https://github.com/user-attachments/assets/08d5f8b0-ad3b-4788-bb27-daa15601ed2a" />
<img width="1394" height="416" alt="image" src="https://github.com/user-attachments/assets/7811e990-c474-49ae-b85c-344e60aafd46" />
5296 potentially good sampling depth?
Samples retained: 27 (OM1), 25 (OM2), 10 (REF)

