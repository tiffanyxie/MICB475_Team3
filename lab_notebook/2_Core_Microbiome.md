# Core Microbiome

## Aim:

To perform core microbiome again at the genus level. To determine what ASVs are unique to each LSTP treatment group (REF/OM1/OM2) and plot those ASVs' relative abundance
Create a Venn diagram using all the ASVs shared and unique to LTSP treatment groups

## Code:

[R Script](../R_scripts/2_CoreMicrobiome.R)

## Results:
*prevelance was set to 0.15 because when it was set to 0.2, it yeilded 0 ASVs for all OM1 and OM2 treatment groups. 

Plot of ASVs, genus level

REF
<img width="1447" height="706" alt="Screenshot 2026-03-04 at 11 33 04 PM" src="https://github.com/user-attachments/assets/b421cb76-96b1-4912-b285-678ecd60c306" />

OM1
<img width="1444" height="704" alt="Screenshot 2026-03-04 at 11 33 28 PM" src="https://github.com/user-attachments/assets/fce54367-0f46-477a-a218-bef9caa900dc" />


OM2
<img width="1441" height="633" alt="Screenshot 2026-03-02 at 12 36 26 AM" src="https://github.com/user-attachments/assets/44a13844-ba95-4b0d-9f04-9b292dbfeecb" />

Venn diagram

When detection was set to 0. Detection set to 0.001 gave the same results.
<img width="636" height="611" alt="Screenshot 2026-03-04 at 11 32 35 PM" src="https://github.com/user-attachments/assets/00976e55-bfe8-42f4-90d0-37c92b8d6951" />

When detection was set to 0.01.

<img width="673" height="671" alt="Screenshot 2026-03-04 at 11 51 56 PM" src="https://github.com/user-attachments/assets/6cf847fd-dae0-471d-8664-eb793541d289" />




*prevalence was set to 0.9 and detection to 0.001.

REF
<img width="1339" height="592" alt="Screenshot 2026-03-11 at 8 52 10 PM" src="https://github.com/user-attachments/assets/ed7cebfd-07bb-4d58-845f-555c2483f2bd" />

OM1
<img width="1443" height="670" alt="Screenshot 2026-03-11 at 8 52 36 PM" src="https://github.com/user-attachments/assets/a4ed2713-75ce-4b2d-a9c3-bbd95e99e1dc" />

OM2
<img width="1442" height="665" alt="Screenshot 2026-03-11 at 8 53 05 PM" src="https://github.com/user-attachments/assets/96b9b66f-7ba1-4cce-b644-7d8c99f13903" />

venn diagram:

<img width="644" height="595" alt="Screenshot 2026-03-11 at 9 26 18 PM" src="https://github.com/user-attachments/assets/ccbb0ccd-0d92-4485-a54d-f0a9c0a90d4e" />

Prevelance set to 0.5 and detection to 0.001
<img width="815" height="615" alt="Screenshot 2026-03-26 at 1 30 19 PM" src="https://github.com/user-attachments/assets/e387d0d0-4907-4c6b-ba51-5f8a3f3eb1be" />

Note that some ASVs were classified into "Incertae sedis"

**Table with the genera used in the venn diagram**

| Reyranella | Incertae_Sedis | Roseiarcus | Acidiferrimicrobium | Mycobacterium | Unknown_Genus (NA)  | Afipia | Edaphobacter | Acidothermus |
| ---------- | -------------- | ---------- | ------------------- | ------------- | --------------------| ------ | ------------ | ------------ |
| OM1        | OM1            | OM1        | OM1                 | OM1           | OM1                 | OM1    |              |              |
| OM2        | OM2            | OM2        | OM2                 | OM2           | OM2                 | OM2    |              |              |
| REF        | REF            | REF        |                     |               | REF                 | REF    | REF          | REF          |

Note that the method used to generate the table combined "Incertae sedis" genera into a single genus, despite there being 4 unique ASVs with differing families or orders. There were also 2 different taxa whose genus is unknown, one from the Xanthobacteraceae family, and one from the Acidobacteriaceae Subgroup 1. 

The R script provides the option of generating a table with the ASV feature id, which yields the 13 unique values reported in the venn diagram.


**Final Thresholds for Manuscript:** 
* 0.5 for prevalence
* 0.001 for abundance