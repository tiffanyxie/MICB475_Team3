# Core Microbiome

## Aim:

To perform core microbiome at the genus level to determine the genera unique to each OM level. To Create a Venn diagram and table to show the genera shared and unique to OM level.

## Code:

[R Script](../R_scripts/2_CoreMicrobiome)

## Results:

Venn diagram 
<img width="661" height="562" alt="Screenshot 2026-04-24 at 2 46 53 PM" src="https://github.com/user-attachments/assets/06435eab-c3cf-4caa-bc5a-14ed484d0cbc" />

Table with unique and shared genera between OM levels

| Genus                 | REF   | OM1   | OM2   |
|----------------------|-------|-------|-------|
| g__Incertae_Sedis    | TRUE  | TRUE  | TRUE  |
| g__Dongia            | TRUE  | TRUE  | TRUE  |
| g__Hypericibacter    | TRUE  | TRUE  | TRUE  |
| g__Robbsia           | TRUE  | TRUE  | TRUE  |
| g__Edaphobacter      | TRUE  | TRUE  | TRUE  |
| g__Reyranella        | TRUE  | TRUE  | TRUE  |
| g__Labrys            | TRUE  | TRUE  | TRUE  |
| g__Roseiarcus        | TRUE  | TRUE  | TRUE  |
| g__Phenylobacterium  | TRUE  | TRUE  | TRUE  |
| g__Acidiferrimicrobium | TRUE | TRUE | TRUE |
| g__Acidothermus      | TRUE  | TRUE  | TRUE  |
| g__Marmoricola       | TRUE  | TRUE  | TRUE  |
| g__Streptomyces      | TRUE  | TRUE  | TRUE  |
| g__Mycobacterium     | TRUE  | TRUE  | TRUE  |
| g__Jatrophihabitans  | TRUE  | TRUE  | TRUE  |
| g__Conexibacter      | TRUE  | TRUE  | TRUE  |
| g__Afipia            | TRUE  | TRUE  | TRUE  |
| g__Inquilinus        | TRUE  | TRUE  | FALSE |
| g__Granulicella      | TRUE  | TRUE  | FALSE |
| g__Baekduia          | TRUE  | FALSE | TRUE  |
| g__Silvibacterium    | TRUE  | FALSE | FALSE |
| g__Rhodopila         | FALSE | TRUE  | TRUE  |
| g__Acidibacter       | FALSE | TRUE  | TRUE  |
| g__Hyphomicrobium    | FALSE | TRUE  | TRUE  |
| g__Bauldia           | FALSE | TRUE  | TRUE  |
| g__Pseudorhodoplanes | FALSE | TRUE  | TRUE  |
| g__Nocardioides      | FALSE | TRUE  | TRUE  |
| g__Solirubrobacter   | FALSE | TRUE  | TRUE  |
| g__Pseudolabrys      | FALSE | TRUE  | TRUE  |
| g__Rhodoplanes       | FALSE | TRUE  | TRUE  |
| g__Variibacter       | FALSE | TRUE  | TRUE  |
| g__Rhodovastum       | FALSE | TRUE  | FALSE |
| g__Pedomicrobium     | FALSE | TRUE  | FALSE |

*TRUE means the genera is present. FALSE means the genera is absent.


**Final Thresholds for Manuscript:** 
* 0.5 for prevalence
* 0.001 for abundance
