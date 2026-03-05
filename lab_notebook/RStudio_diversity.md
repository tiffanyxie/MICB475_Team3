# Initial Diversity analyses

## Aim:
* Support experimental aims and assess important variables
* Generate phyloseq object and rarefied phyloseq object
* Perform alpha and beta diversity statistical analyses

## Code:
[R Script](soil_export/Diversity_by_metadata.R)

## Alpha Diversity Results:
### OM Removal
**IDFBC Ecozone**

<img width="400" height="400" alt="a1e7f66e-5d57-41ef-a3f6-813d74a699f5" src="https://github.com/user-attachments/assets/9498df5c-e916-4149-8010-5262d15cf5d2" />
<img width="400" height="400" alt="2a0a0217-99ef-4fa0-a8b4-c1a4c7e73050" src="https://github.com/user-attachments/assets/77cf663c-3a1b-40a2-bde6-f3bbff7ff106" />

No Significant Differences

**SBSBC Ecozone**

<img width="400" height="400" alt="080c8412-4a2e-412c-8282-3e77b0dd1b34" src="https://github.com/user-attachments/assets/b57c596a-1ac4-4976-9464-c8e7c80ea2aa" />
<img width="400" height="400" alt="00f3b995-3e35-46ad-bc84-927560dcf43d" src="https://github.com/user-attachments/assets/51d9e67b-ce04-49a3-9913-c467522d3671" />

No Significant Differences

**All of BC**

<img width="400" height="400" alt="1f1e35b6-f421-42f5-a005-fa34a57190b7" src="https://github.com/user-attachments/assets/ae93598e-f242-4024-91a6-386868217cde" />
<img width="400" height="400" alt="b7346ecc-64e3-48b6-8e1d-a256da81058b" src="https://github.com/user-attachments/assets/50a2d880-b12d-46ea-9ea3-4c632bf09200" />

No Significant Differences

### Compaction Treatment
**IDFBC Ecozone**

<img width="400" height="400" alt="ec33ee0f-bfea-4813-832c-d30574036f0e" src="https://github.com/user-attachments/assets/8f8cf61d-61f3-4faf-a66f-88a26372ffa9" />
<img width="400" height="400" alt="1bde320f-09f7-4a10-bd62-0cef61fe35c2" src="https://github.com/user-attachments/assets/0bb4b372-ea49-4ce5-9ed6-bc1964cb70a9" />

Significance: 
* Chao1 - p = 0.04843, not identified by ANOVA, likely REF vs C0
* ACE - p = 0.04276, not identified by ANOVA, likely REF vs C0

**SBSBC Ecozone**

<img width="400" height="400" alt="cef1573b-cb1f-4869-9c4b-e5fe49de5673" src="https://github.com/user-attachments/assets/ad47443d-1c28-40e0-810f-8fd578e21331" />
<img width="400" height="400" alt="f7bebb96-fc0e-4da0-8c09-bcd7bc45d018" src="https://github.com/user-attachments/assets/cee2d418-dd4c-499b-bddb-b0b6dd5ed75d" />

No Significant Differences

**All of BC**

<img width="400" height="400" alt="9b09fc43-dc7a-484e-b8f8-a1f6258fe916" src="https://github.com/user-attachments/assets/cc4cd8da-719c-40a2-a480-872f9efc04df" />
<img width="400" height="400" alt="0055ef03-3436-45a3-9691-07e956fbc4fb" src="https://github.com/user-attachments/assets/f2dc51fb-f406-4daf-a025-58ab3cd8884b" />

No Significant Differences (Shocker, right?)

### Site
**IDFBC Ecozone**

<img width="400" height="400" alt="330afd9f-50ec-4bf1-baaf-e1046757016f" src="https://github.com/user-attachments/assets/1c382f6b-3360-46d6-9207-271c8e5bbc2d" />
<img width="400" height="400" alt="fd71ab83-4e4a-44c3-945a-c36d8768bcbc" src="https://github.com/user-attachments/assets/1bad6ba7-0ddd-4d7e-bf0b-1a132389aa1b" />

Significance:
* PD - p = 0.03009, Comparison of OC to DC yields 0.046 from ANOVA

**SBSBC Ecozone**

<img width="400" height="400" alt="97a4d234-62c7-4f99-8fef-005d65810e52" src="https://github.com/user-attachments/assets/5012cefb-eeff-4940-b093-e365b7b4883c" />
<img width="400" height="400" alt="ddaf717c-d2e4-4549-b4c8-adfd5a1ded90" src="https://github.com/user-attachments/assets/8aed1c0c-be71-4786-a882-dc3a3f3a8116" />

No Significant Differences

**All BC**

<img width="400" height="400" alt="18d27789-bc85-4a02-b9c5-046534f39ff5" src="https://github.com/user-attachments/assets/c75ba7a3-686d-46a3-ba4e-f2e6eda13de6" />
<img width="400" height="400" alt="8d877d4f-b3c9-4a26-82a7-a9018f594938" src="https://github.com/user-attachments/assets/148b3a55-20ee-4745-8c84-d2f49121510e" />

Significant differences across all metrics, as expected

### Moisture Content

**IDFBC**

<img width="400" height="400" alt="9c4ffee9-ab7c-41b8-af7a-e42977e18366" src="https://github.com/user-attachments/assets/da7a22ff-af61-4271-838a-2dd960db8d8b" />

No Significant Correlation

**SBSBC**

<img width="400" height="400" alt="b1719b1d-5fff-4e34-801b-6ceabe77084e" src="https://github.com/user-attachments/assets/794ee182-7c74-43ad-aa5c-8030456f7784" />

No Significant Correlation

**All BC**

<img width="400" height="400" alt="94f77f38-742f-427b-acc3-cc9dc84fbbd6" src="https://github.com/user-attachments/assets/48bf6a0b-c83f-4523-b54b-0174b95e6c98" />

Significant Correlation across all metrics, Moisture across ecozones may vary significantly

### Total Carbon

**IDFBC**

<img width="400" height="400" alt="42e11382-f6ee-4117-a030-382d6710509a" src="https://github.com/user-attachments/assets/6c96c410-895e-4e67-80fe-ca7ca46da2f4" />

No Significant Correlation

**SBSBC**

<img width="400" height="400" alt="410b44a1-ba33-4ddd-b4e8-d3eaf2aa267c" src="https://github.com/user-attachments/assets/f919f168-e1b1-44ad-81ce-a6d726fbb0ef" />

No Significant Correlation

**All BC**

<img width="400" height="400" alt="9929a7ec-2d56-4d85-91ed-c64f24f72083" src="https://github.com/user-attachments/assets/b55e3832-54f5-4981-b3b3-e5948a5b931e" />

Significant Correlation across all metrics, Total Carbon across ecozones may vary significantly

### Total Nitrogen

**IDFBC**



**SBSBC**



**All BC**



### CN Ratio

**IDFBC**



**SBSBC**



**All BC**


### pH

**IDFBC**



**SBSBC**



**All BC**


### Soil Bulk Density

**IDFBC**



**SBSBC**



**All BC**





