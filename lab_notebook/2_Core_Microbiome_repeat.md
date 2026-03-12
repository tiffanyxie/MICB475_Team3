# Core Microbiome

## Aim:

To perform core microbiome again at the genus level. To determine what ASVs are unique to each LSTP treatment group (REF/OM1/OM2) and plot those ASVs' relative abundance
Create a Venn diagram using all the ASVs shared and unique to LTSP treatment groups

## Code:

[R Script](https://github.com/tiffanyxie/MICB475_Team3/blob/55c54bfa1f9156a4a6b3ba350240fb495d8dd341/R_scripts/2_coremicrobiome_repeat.R)

## Results:

*prevalence was set to 0.9 and detection to 0.001.

REF
<img width="1339" height="592" alt="Screenshot 2026-03-11 at 8 52 10 PM" src="https://github.com/user-attachments/assets/ed7cebfd-07bb-4d58-845f-555c2483f2bd" />

OM1
<img width="1443" height="670" alt="Screenshot 2026-03-11 at 8 52 36 PM" src="https://github.com/user-attachments/assets/a4ed2713-75ce-4b2d-a9c3-bbd95e99e1dc" />

OM2
<img width="1442" height="665" alt="Screenshot 2026-03-11 at 8 53 05 PM" src="https://github.com/user-attachments/assets/96b9b66f-7ba1-4cce-b644-7d8c99f13903" />

venn diagram:

<img width="644" height="595" alt="Screenshot 2026-03-11 at 9 26 18 PM" src="https://github.com/user-attachments/assets/ccbb0ccd-0d92-4485-a54d-f0a9c0a90d4e" />

Note that some ASVs were classified into "Incertae sedis"

**Table with the genera used in the venn diagram**

| Reyranella | Incertae_Sedis | Roseiarcus | Acidiferrimicrobium | Mycobacterium | Unknown_Genus (NA)  | Afipia | Edaphobacter | Acidothermus |
| ---------- | -------------- | ---------- | ------------------- | ------------- | --------------------| ------ | ------------ | ------------ |
| OM1        | OM1            | OM1        | OM1                 | OM1           | OM1                 | OM1    |              |              |
| OM2        | OM2            | OM2        | OM2                 | OM2           | OM2                 | OM2    |              |              |
| REF        | REF            | REF        |                     |               | REF                 | REF    | REF          | REF          |

Note that the method used to generate the table combined "Incertae sedis" genera into a single genus, despite there being 4 unique ASVs with differing families or orders. There were also 2 different taxa whose genus is unknown, one from the Xanthobacteraceae family, and one from the Acidobacteriaceae Subgroup 1. 

The R script provides the option of generating a table with the ASV feature id, which yields the 13 unique values reported in the venn diagram.
