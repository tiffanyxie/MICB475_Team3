# Compare soil conditions between OM Removal

## Aim
* Compare soil conditions (moisture, total carbon, total nitrogen, C/N ratio, and pH)  between OM removal

## Code

[R Script](../R_scripts/1_Metadata_Analysis.R)

## Results
**Overview**
* Moisture content lower in OM1 and OM2 vs REF in SBSBC
* Total carbon decreases with increasing OM removal in overall BC and IDFBC, but no significance in SBSBC
* Lower nitrogen in OM removal samples in BC, no significance when subset to ecozones
* Lower CN ratio in OM2 vs OM1 and REF in BC and IDFBC, no significance in SDSBC
* pH higher in OM2 vs REF in IDFBC only

### Moisture Content

<img width="898" height="309" alt="image" src="https://github.com/user-attachments/assets/0d5d83de-6043-4c8c-9102-8e8cf49bf402" />

**BC**
* No significant difference across OM treatments within BC
Kruskal-Wallis: chi-squared = 5.2861, df = 2, p-value =
0.07114

Dunn test

| Comparison | Z         | P.unadj    | P.adj     |
|------------|-----------|------------|-----------|
| OM1 - OM2  | 0.2991179 | 0.76485007 | 0.76485007|
| OM1 - REF  | -2.0285439| 0.04250477 | 0.08500954|
| OM2 - REF  | -2.2357807| 0.02536614 | 0.07609842|


**IDFBC**
* No significant difference
Kruskal-Wallis chi-squared = 0.55683, df = 2, p-value = 0.757
  
Dunn test

| Comparison   | Z          | P.unadj   | P.adj     |
|--------------|------------|-----------|-----------|
| OM1 - OM2    | -0.4012725 | 0.6882195 | 0.6882195 |
| OM1 - REF    | -0.7303742 | 0.4651615 | 1.0000000 |
| OM2 - REF    | -0.4466317 | 0.6551410 | 1.0000000 |

**SBSBC**

* OM1 and OM2 significantly lower moisture than REF
Kruskal-Wallis chi-squared = 8.889, df = 2, p-value = 0.01174
  
Dunn test

| Comparison   | Z          | P.unadj     | P.adj      |
|--------------|------------|-------------|-------------|
| OM1 - OM2    |  0.6496412 | 0.515924036 | 0.515924036 |
| OM1 - REF    | -2.4980307 | 0.012488538 | 0.024977076 |
| OM2 - REF    | -2.9476548 | 0.003201944 | 0.009605832 |

### Total Carbon
<img width="716" height="309" alt="image" src="https://github.com/user-attachments/assets/53a9bead-4d1b-4827-9e58-fa236f23381f" />

**BC**
* Lower carbon with increasing severity of OM removal
Kruskal-Wallis chi-squared = 17.252, df = 2, p-value = 0.0001794

Dunn test
  
| Comparison |       Z   |    P.unadj    |    P.adj      |
|------------|-----------|---------------|---------------|
| OM1 - OM2  |  2.077556 | 3.775026e-02  | 0.0377502559  |
| OM1 - REF  | -2.636675 | 8.372303e-03  | 0.0167446064  |
| OM2 - REF  | -4.102414 | 4.088615e-05  | 0.0001226585  |

**IDFBC**
* Lower carbon with increasing severity of OM removal
Kruskal-Wallis chi-squared = 20.093, df = 2, p-value = 4.334e-05

Dunn test

| Comparison |       Z   |    P.unadj    |    P.adj      |
|------------|-----------|---------------|---------------|
| OM1 - OM2  |  2.210750 | 2.705313e-02  | 2.705313e-02  |
| OM1 - REF  | -2.865934 | 4.157813e-03  | 8.315625e-03  |
| OM2 - REF  | -4.429170 | 9.459630e-06  | 2.837889e-05  |


**SBSBC**
* Not significant

Kruskal-Wallis chi-squared = 5.1755, df = 2, p-value = 0.07519

  Dunn test
  

| Comparison |       Z   |   P.unadj   |   P.adj    |
|------------|-----------|-------------|------------|
| OM1 - OM2  |  0.9413612 | 0.3465198  | 0.3465198 |
| OM1 - REF  | -1.6102426 | 0.1073449  | 0.2146898 |
| OM2 - REF  | -2.2713616 | 0.0231251  | 0.0693753 |

### Total Nitrogen
<img width="751" height="309" alt="image" src="https://github.com/user-attachments/assets/dad8d693-14dc-4e96-9db6-bd7c535094be" />

**BC**
* Lower nitrogen in OM2 vs ref

Kruskal-Wallis chi-squared = 6.9956, df = 2, p-value =
0.03026
  
Dunn test

| Comparison |       Z    |   P.unadj    |   P.adj     |
|------------|------------|--------------|-------------|
| OM1 - OM2  |  0.5859777 | 0.557890484  | 0.55789048  |
| OM1 - REF  | -2.2085421 | 0.027206506  | 0.05441301  |
| OM2 - REF  | -2.6178121 | 0.008849552  | 0.02654866  |

**IDFBC**
* No significance

Kruskal-Wallis chi-squared = 5.4791, df = 2, p-value = 0.0646

Dunn test
| Comparison |       Z    |   P.unadj    |   P.adj     |
|------------|------------|--------------|-------------|
| OM1 - OM2  |  0.3742875 | 0.70819045   | 0.70819045  |
| OM1 - REF  | -2.0290693 | 0.04245123   | 0.08490246  |
| OM2 - REF  | -2.2937305 | 0.02180598   | 0.06541795  |

**SBSBC**
* No significance
  
Kruskal-Wallis chi-squared = 3.791, df = 2, p-value = 0.1502

Dunn test
| Comparison |       Z    |   P.unadj    |   P.adj     |
|------------|------------|--------------|-------------|
| OM1 - OM2  |  0.2818725 | 0.77804129   | 0.7780413   |
| OM1 - REF  | -1.7058292 | 0.08803988   | 0.1760798   |
| OM2 - REF  | -1.8979476 | 0.05770299   | 0.1731090   |

### CN Ratio
<img width="796" height="309" alt="image" src="https://github.com/user-attachments/assets/8f8051cc-84c1-47dc-bf32-d12174c1d31b" />

**BC**
* Significantly lower CN ratio in OM2 vs OM1 and REF
  
Kruskal-Wallis chi-squared = 14.1836, df = 2, p-value = 0.0008319
Dunn test
| Comparison |       Z     |    P.unadj      |    P.adj       |
|------------|------------|----------------|----------------|
| OM1 - OM2  |  2.727485  | 0.0063819073   | 0.012763815    |
| OM1 - REF  | -1.472342  | 0.1409284542   | 0.140928454    |
| OM2 - REF  | -3.402241  | 0.0006683563   | 0.002005069    |

**IDFBC**
* Significantly lower CN ratio in OM2 vs OM1 and REF

Kruskal-Wallis chi-squared = 8.3385, df = 2, p-value = 0.01546

Dunn test
| Comparison |       Z     |   P.unadj    |   P.adj     |
|------------|------------|--------------|-------------|
| OM1 - OM2  |  2.2696563 | 0.02322844   | 0.04645688  |
| OM1 - REF  | -0.8675078 | 0.38566385   | 0.38566385  |
| OM2 - REF  | -2.4723972 | 0.01342103   | 0.04026309  |

**SBSBC**
* No significant difference
Kruskal-Wallis chi-squared = 4.5462, df = 2, p-value = 0.103

Dunn test

| Comparison |       Z    |   P.unadj    |   P.adj     |
|------------|------------|--------------|-------------|
| OM1 - OM2  |  1.384675  | 0.16615202   | 0.3323040   |
| OM1 - REF  | -1.033792  | 0.30123353   | 0.3012335   |
| OM2 - REF  | -2.012621  | 0.04415449   | 0.1324635   |

### pH
<img width="860" height="438" alt="image" src="https://github.com/user-attachments/assets/d34173e2-a56a-461f-a09c-b1e7ba1ad857" />


**BC**
* No significance
Kruskal-Wallis chi-squared = 3.3378, df = 2, p-value = 0.1885

Dunn test
| Comparison | Z | P.unadj | P.adj |
|---|---|---|---|
| OM1 - OM2 | -0.03427439 | 0.97265834 | 0.9726583 |
| OM1 - REF | 1.69798354 | 0.08951087 | 0.1790217 |
| OM2 - REF | 1.71748420 | 0.08589074 | 0.2576722 |

**IDFBC**
* Significantly higher in OM2 vs REF
Kruskal-Wallis chi-squared = 7.3333, df = 2, p-value = 0.02556

Dunn test

| Comparison | Z | P.unadj | P.adj |
|---|---|---|---|
| OM1 - OM2 | -0.7385489 | 0.460180935 | 0.46018094 |
| OM1 - REF | 2.1759707 | 0.029557452 | 0.05911490 |
| OM2 - REF | 2.6982037 | 0.006971477 | 0.02091443 | 

**SBSBC**
* No significance
Kruskal-Wallis chi-squared = 2.6324, df = 2, p-value = 0.2681

Dunn test
| Comparison |       Z       |   P.unadj    |   P.adj     |
|------------|---------------|--------------|-------------|
| OM1 - OM2  |  0.5883648    | 0.5562875    | 0.5562875   |
| OM1 - REF  |  1.6223313    | 0.1047324    | 0.3141973   |
| OM2 - REF  |  1.1965743    | 0.2314725    | 0.4629451   |


