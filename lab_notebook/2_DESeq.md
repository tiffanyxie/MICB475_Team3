# DESeq

## Aim:
* Determine ASVs (or Taxa) that have decreased or increased in presence relative to a reference

## Code:
[R Script](../R_scripts/2_DESeq.R)

## Results:
Volcano Plots
<img width="1124" height="884" alt="image" src="https://github.com/user-attachments/assets/d2813a7f-48fb-4879-b9aa-6cd87efddc69" />
<img width="1124" height="884" alt="image" src="https://github.com/user-attachments/assets/8baef78f-f5e5-4b92-9e7b-fe6d8c3a7f74" />
<img width="1124" height="884" alt="image" src="https://github.com/user-attachments/assets/db470b8d-690a-494f-8a37-954a809ebceb" />

Taxa bar plots
<img width="1370" height="884" alt="image" src="https://github.com/user-attachments/assets/18c42d71-e939-46ad-beb0-12abde0193cd" />
<img width="1370" height="884" alt="image" src="https://github.com/user-attachments/assets/fce3fe13-f39a-4817-96fd-1721ce43714a" />
<img width="1370" height="884" alt="image" src="https://github.com/user-attachments/assets/ed96ec13-4a5a-4504-8d39-3503e74004be" />

### Re-do with glom | March 9, 2026
Had to adjust padj<0.05 & abs(log2FoldChange)>1 from padj<0.01 & abs(log2FoldChange)>2, or else nothing shows
<img width="954" height="882" alt="image" src="https://github.com/user-attachments/assets/03f46a61-700f-4e06-acbf-4fc6f79f99f8" />
<img width="954" height="882" alt="image" src="https://github.com/user-attachments/assets/47b74188-0c81-4ef4-922c-9a12488048a9" />
<img width="954" height="882" alt="image" src="https://github.com/user-attachments/assets/91a4ba6f-f284-4937-a7eb-ca0d2e87a96e" />

<img width="1172" height="882" alt="image" src="https://github.com/user-attachments/assets/43ed83f5-4afc-4178-8a8b-6e6033c7206d" />
<img width="1172" height="882" alt="image" src="https://github.com/user-attachments/assets/6952e16e-902d-4e39-8ba0-370fd3c3733f" />
<img width="1172" height="882" alt="image" src="https://github.com/user-attachments/assets/990a902f-a436-4c61-a5b7-887ca3dd9bb1" />
