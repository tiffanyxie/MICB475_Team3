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

<img width="400" height="400" alt="c9f0cf94-1033-4b29-a0b1-5c40780f7db2" src="https://github.com/user-attachments/assets/a2310197-b1c4-462e-bd39-13f65a49c9d8" />

No Significant Correlation

**SBSBC**

<img width="400" height="400" alt="51ea213e-c114-4465-b337-cd2bb21f0255" src="https://github.com/user-attachments/assets/33f9f9ef-9014-4e8c-9441-c54a5c0667ad" />

No Significant Correlation

**All BC**

<img width="400" height="400" alt="c6bcaa6c-a0c3-4a63-94d6-c00b20b920c2" src="https://github.com/user-attachments/assets/c160f38a-73f2-4fb9-8319-c47a217711c6" />

Significant Correlation across all metrics, Total Nitrogen across ecozones may vary significantly

### CN Ratio

**IDFBC**

<img width="400" height="400" alt="f713c0d2-08d2-43dd-b6a8-04a91e7ce6f2" src="https://github.com/user-attachments/assets/6756af8c-e50d-4fc5-9e18-0f1d01ab24d0" />

No Significant Correlation

**SBSBC**

<img width="400" height="400" alt="683ebd9c-c707-4073-aca8-101702bf6821" src="https://github.com/user-attachments/assets/8c1f04e3-1834-4cc0-b91c-53d831dea5f2" />

No Significant Correlation

**All BC**

<img width="400" height="400" alt="bf77c711-3ff8-4f22-9c78-b39ff437ddc5" src="https://github.com/user-attachments/assets/5d41ef67-a34f-4233-b267-3ccc9c4dd788" />

No Significant Correlation

### pH

**IDFBC**

<img width="400" height="400" alt="3fff1774-fa77-49fd-8dd3-3fcfceb094b3" src="https://github.com/user-attachments/assets/31a656f5-0051-4068-ad6c-d9993916b092" />

No Significant Correlation

**SBSBC**

<img width="400" height="400" alt="24179bb1-bc8d-44cf-b673-b296ac8b56aa" src="https://github.com/user-attachments/assets/40359433-7ad2-4124-ad13-62b05fbb5f5f" />

Significance
* PD - p = 0.006103

**All BC**

<img width="400" height="400" alt="17aa6efe-6050-4e75-ba64-429dbdfb72e5" src="https://github.com/user-attachments/assets/985f8ff7-77b9-47ca-bb28-16454b230a5f" />
<img width="400" height="400" alt="050dd5e0-92a6-4094-b2a2-9c941667e4c5" src="https://github.com/user-attachments/assets/dd22131a-14e5-4804-83a4-c07f1e2f3ebc" />

Significant Correlation across all metrics, pH across ecozones may vary significantly. Box plot shows the difference in pH across Ecozones

### Soil Bulk Density

**IDFBC**

<img width="400" height="400" alt="6b9d58e7-7af1-4ed7-9c01-3c554ba817f3" src="https://github.com/user-attachments/assets/a780055f-7205-4cbb-a2d5-4c9a2835f8f1" />

No Significant Correlation

**SBSBC**

<img width="400" height="400" alt="54dfd76b-2d4c-4a0e-ac56-626b57d53caf" src="https://github.com/user-attachments/assets/b18eeb34-9050-445e-a825-f9836b83bf04" />

No Significant Correlation

**All BC**

<img width="400" height="400" alt="19d3314c-ca5b-461d-aa25-0b55b49661c4" src="https://github.com/user-attachments/assets/5f6c36f7-8a6c-4bf1-a3b8-88303efe7ac5" />

Significant Correlation across all metrics, soil bulk density across ecozones may vary significantly.

## Beta Diversity Results

**IDFBC**

<img width="400" height="400" alt="a7a1d9cf-44fb-4c27-8084-e40f5857d0c0" src="https://github.com/user-attachments/assets/bdb6067e-428e-4d43-ae9a-cbf0036a5f37" />
<img width="400" height="400" alt="cf656d6f-7ead-47fa-8cee-a979daac928c" src="https://github.com/user-attachments/assets/86ea0b1a-4ae2-4d81-a0f6-ec5ba16d1711" />


**SBSBC**

<img width="400" height="400" alt="bdd69c08-2a90-4721-b3dc-c3bee13a854a" src="https://github.com/user-attachments/assets/dabcd7b4-910f-48f0-babe-51c5567ea23c" />
<img width="400" height="400" alt="d65a91c7-7b4b-407c-8699-e4662df0ac17" src="https://github.com/user-attachments/assets/70291236-002a-461f-9379-d18f702225a3" />



**All BC**

<img width="400" height="400" alt="ea397fc0-262d-4dc8-9f93-0512f02aa183" src="https://github.com/user-attachments/assets/d13ecbfb-f30b-49bd-8490-c99897ab0d83" />
<img width="400" height="400" alt="4bc5ad9b-f773-4b0a-a9e9-6dbf8fbb5c89" src="https://github.com/user-attachments/assets/966af4c3-8614-4b8f-a065-344b901eca83" />


