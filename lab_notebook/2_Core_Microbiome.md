# Core Microbiome

## Aim:
* Determine what ASVs are unique to each LSTP treatment group (REF/OM1/OM2) and plot those ASVs' relative abundance
* Create a Venn diagram using all the ASVs shared and unique to LTSP treatment groups 

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

Prevelance set to 0.5 and detection to 0.001
<img width="815" height="615" alt="Screenshot 2026-03-26 at 1 30 19 PM" src="https://github.com/user-attachments/assets/e387d0d0-4907-4c6b-ba51-5f8a3f3eb1be" />




