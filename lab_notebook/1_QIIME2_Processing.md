# 1. QIIME2 Data Processing

## Aim:
* Import and demultiplex data from soil dataset
* Denoise data by trimming and removing low quality reads
* Cluster unique reads into amplicon sequence variants (ASVs)
* Filter out mitochondrial and chloroplast DNA

## Code
[QIIME2 Script](../QIIME2/qiime2_script.sh)

## Results

### Post-demultiplexing
[demux.qzv](../QIIME2/output/demux.qzv)

**Demux Sequence Statistics**:

![Sequence counts summary](../QIIME2/output/demux_seq_counts.png)

**Quality Score per Base**:

![Base quality score](../QIIME2/output/demux_base_qual.png)

**rep-seqs.qzv**
<img width="627" height="197" alt="image" src="https://github.com/user-attachments/assets/443c6ff3-e1f5-44b3-8a23-9d4d4109e96e" />


**Demux Sequence Lengths**:

![Sequence lengths](../QIIME2/output/demux_seq_length.png)

Decided to trim to 390 bp

### ASVs



<img width="870" height="590" alt="image" src="https://github.com/user-attachments/assets/744cd9fc-a4db-405b-87fc-26fec3fe3494" />

### Taxonomy Analysis
* Trained classifier using universal 27F and 519R primers covering V1 to V3
* Performed taxonomic analysis using this classifier
![Taxa bar plot](../QIIME2/output/level-4-bars.svg)

### Mitochondria and chloroplast filtered table
* Lost 1 sample and ~ 7000 ASVs
<img width="699" height="555" alt="image" src="https://github.com/user-attachments/assets/30485b64-96a9-46a3-a114-fbb8b75d5353" />

### Organic layer + no herbicide + BC only filtering
<img width="709" height="566" alt="image" src="https://github.com/user-attachments/assets/08d5f8b0-ad3b-4788-bb27-daa15601ed2a" />
<img width="1394" height="416" alt="image" src="https://github.com/user-attachments/assets/7811e990-c474-49ae-b85c-344e60aafd46" />
5296 potentially good sampling depth?
Samples retained: 27 (OM1), 25 (OM2), 10 (REF)

### Rarefaction curve

<img width="1088" height="462" alt="Screenshot 2026-02-18 at 12 56 29 PM" src="https://github.com/user-attachments/assets/e0fd7861-8b08-46ae-bcd2-61bb493bfcc0" />
Sequencing depth of ~3000 identified as the minimum

<img width="1431" height="422" alt="Screenshot 2026-02-18 at 12 58 10 PM" src="https://github.com/user-attachments/assets/60dec78d-f1d5-466c-a55c-3d672b1a1a81" />
Sampling depth revision to 3310 <br\>

Samples retained: 45 (OM1), 41 (OM2), 15 (REF)

* Confirmed with Bessie, will be used going forward
 
