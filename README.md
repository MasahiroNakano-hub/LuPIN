# LuPIN (Lupus-derived Peripheral immune cell Integrated panels)

## Reference: 
### Nakano M et al. Deep immune profiling pinpoints cellular and molecular drivers of lupus immunopathology (Preprint: https://www.biorxiv.org/content/10.64898/2025.12.08.693110v2)


## Overview:
### ・To pinpoint new disease-relevant cell states and their molecular profiles, we performed an in-depth investigation of multimodal single-cell datasets comprising ~2.1 million peripheral blood mononuclear cells from 346 donors.
### ・By resolving 123 fine-grained cell states across 27 cell types, we identified previously uncharacterized populations distinctively associated with clinical severity and treatment status, including GZMK+GZMH+HLA-DR+ effector memory CD8+ T cells (double-positive [DP] EMCD8) and FOXO1+ARHGAP15+ T cells.
### ・Through extensive statistical frameworks and multimodal approaches, we delineated their aberrant immune signaling networks, transcriptional regulators, key surface proteins, T cell receptor repertoires, and genetic/epigenetic landscapes, underscoring them as candidate drivers of SLE immunopathology.

<img width="1644" height="867" alt="スクリーンショット 2026-03-11 14 04 58" src="https://github.com/user-attachments/assets/95f8b5e7-c3a7-4767-b07a-5010b4283690" />


## LuPIN Resource
### 1: High-resolution immune cell reference panel
#### lineage-level & cell-type-level Symphony reference objects
#### Researchers can map their datasets onto LuPIN for precise cell-state annotation and mechanistic inference.

### 2: A comprehensive atlas of ligand-receptor signaling networks across 123 cell states
####  ligand-receptor pairs with "sca.LRscore" > 0.6 from LIANA.

### 3: Multimodal feature atlas across fine-grained cell states.
#### Cell-state-specific surface markers within each cell type.
#### Cell-state-specific TCR V genes/clonotypes.


## Results:
### 1: Definition of 123 robust cell states using fine-grained clustering. (Fig.2)
#### To identify granular cell states of SLE, we employed a three-step intensive clustering approach for PBMCs (lineage; cell types; cell states). Louvain clustering on 27 cell types resolved 123 fine-grained cell states with high reproducibility in independent datasets.
<img width="810" height="912" alt="スクリーンショット 2026-03-16 23 11 00" src="https://github.com/user-attachments/assets/b1626c11-04b6-4d71-9dfc-e53f0f502b34" />

### 2: Mapping cell state compositional changes in SLE. (Fig.3)
#### We assessed cell-state compositional changes within each cell type using CNA. 29 cell states, including known (e.g., Tph and ABCs) and previously uncharacterized ones (e.g., GZMK+GZMH+HLA-DRB1+ [DP] EMCD8 and FOXO1+ARHGAP15+ T cells), showed evidence of remarkable expansion especially in severe SLE.
<img width="808" height="810" alt="スクリーンショット 2026-03-16 23 12 36" src="https://github.com/user-attachments/assets/bfcb54c2-83ee-4468-9e16-8cd9c200abe1" />

### 3: Statistical framework to dissect quantitative and qualitative changes within cell types. (Fig.4)
#### To address the limitation of cell-type-level analysis, we constructed a generalized linear mixed-effect model jointly regressing each single-cell expression by (i) a phenotype-driven cell-state abundance term (i.e., CNA neighborhood correlation) and (ii) a direct phenotype term capturing per-cell gene dysregulation. This model demonstrated that key immunological pathways and SLE risk variants were predominantly enriched in cell-state compositional changes within cell types.
<img width="809" height="784" alt="スクリーンショット 2026-03-16 23 21 11" src="https://github.com/user-attachments/assets/8cee9629-9c74-4bb0-aca8-d05320b4b135" />

### 4: Characterization of aberrant immune networks between disease-relevant cell states. (Fig.5)
#### Hierarchical clustering based on cell-state proportions revealed five distinct sample groups (SGs) and six cell-state groups (CGs). CG2 included both established pathogenic and newly identified cell populations. SG1-2 with remarkable expansion of CG2 comprised a higher proportion of severe phenotypes. Focusing on CG2, ligand-receptor network analysis suggested that these newly identified cell populations engage in intricate immune crosstalk with established pathogenic states.
<img width="810" height="616" alt="スクリーンショット 2026-03-16 23 33 18" src="https://github.com/user-attachments/assets/d5014b0f-be99-42a7-93ff-c904d8d16494" />

### 5: FOXO1 is a candidate driver of the pathogenic CD4+ T cell state gene program. (Fig.6)
#### We characterized the multimodal molecular profiles (GWAS candidate genes/epigenome/surface markers/TCR repertoires) in newly identified CD4+ T cell states within CG2. Specifically, FOXO1+ARHGAP15+ cells showed the multiple evidence of TCR signaling activation. Both computational and experimental approaches identified FOXO1 as a transcriptional regulator candidate of this cell state, which could be a new therapeutic target.
<img width="808" height="696" alt="スクリーンショット 2026-03-16 23 33 43" src="https://github.com/user-attachments/assets/113b0bd5-6bbb-4d9b-818a-acd71c3c1b56" />

### 6: DP EMCD8 cells as a candidate pathogenic player in SLE. (Fig.7)
#### In-depth investigation pinpointed DP EMCD8 (CG2) as the cellular origin of multiple pro-inflammatory cytokines (IFNG, CD70, CCL5, IL16). Specifically, DP EMCD8 and Tph shared several transcriptional signatures and they may synergistically drive the maturation of pathogenic B cells (ABCs). DP EMCD8 showed the evidence of TCR clonal expansion. Furhermore, multimodal profiling indentified four key surface markers of this population, which could be a new therapeutic target.
<img width="808" height="642" alt="スクリーンショット 2026-03-16 23 34 03" src="https://github.com/user-attachments/assets/a21873d7-6bdc-4b26-adb5-1aa3f36d8acd" />


## Systems, softwares, and installation
#### R softwares: Seurat (v.5.1.0), harmony (v.1.2.3), symphony (v.0.1.1), lme4 (v.1.1-36), clusterProfiler (v.4.14.0), liana (v.0.1.10), decoupleR (v.2.16.0), monocle3 (v.1.4.26) etc.
#### python softwares: cna (v.0.1.6), ldsc (v1.0.1), etc.
#### Environments can be prepared by Anaconda (~30 min).　
#### See each code and demo for more details.
