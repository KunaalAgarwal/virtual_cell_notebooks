# Data Access Details for Public Perturbation Datasets

All datasets below use **human** cells. Every cell line is a human-derived immortalized line or primary cell culture.

## Current Project Infrastructure

The project uses **scanpy/anndata** (h5ad format) and streams large files via **h5py**. Installed: `scanpy`, `gcsfs`, `pandas`, `scipy`. **Not installed**: `pertpy`, `huggingface_hub`, `torch`, `GEARS`. GCS access is anonymous (public buckets only). No BigQuery setup.

## Cell Line Reference

| Cell Line | Organism | Tissue of Origin | Disease Context | Notes |
|-----------|----------|-------------------|-----------------|-------|
| K562 | Human | Bone marrow | Chronic myelogenous leukemia = blood cancer (CML) | Derived from a 53-year-old female CML patient in blast (progenitor of blood cells) crisis. One of the most widely used cell lines in biology. Suspension culture (non-adherent). Has the BCR-ABL (couple of genes are fused which creates the mutation that leads to the cancer) |
| RPE1 (hTERT-RPE1) | Human | Retina | Non-cancerous (immortalized normal) | Retinal pigmented epithelial cells immortalized with telomerase (hTERT). Near-diploid, adherent. One of the few widely used non-cancer human cell lines. Originally from a healthy donor. |
| HepG2 | Human | Liver | Hepatocellular carcinoma (liver cancer) | Derived from a 15-year-old male with well-differentiated hepatocellular carcinoma. Adherent, epithelial morphology. Widely used for drug metabolism and hepatotoxicity studies because it retains many liver-specific functions (albumin secretion, drug-metabolizing enzymes). |
| Jurkat | Human | Blood (T-cells) | Acute T-cell leukemia | Derived from a 14-year-old male with T-cell acute lymphoblastic leukemia. Suspension culture. The standard model for T-cell receptor signaling and immune cell biology. |
| iPSC lines (34 lines) | Human | Reprogrammed (various donors) | Non-cancerous (reprogrammed pluripotent) | Induced pluripotent stem cells from 26 healthy donors across diverse genetic backgrounds. Pluripotent — can differentiate into any cell type. The closest available model to the H1 hESC line used in VCC 2025. |
| A549 | Human | Lung | Non-small cell lung adenocarcinoma (specific kind of lung cancer)| Derived from a 58-year-old male. Adherent, epithelial. Widely used in lung cancer, drug screening, and respiratory virus research. Expresses many lung-specific markers. |
| MCF7 | Human | Breast | Breast adenocarcinoma (specific kind of breast cancer) | Derived from a 69-year-old female. Adherent, estrogen receptor-positive (ER+). The most-studied breast cancer cell line. |
| HT29 | Human | Colon | Colorectal adenocarcinoma | Derived from a 44-year-old female. Adherent, epithelial. Can differentiate into goblet-like or absorptive enterocyte-like cells under certain conditions. Standard colorectal cancer model. |
| HAP1 | Human | Bone marrow | Derived from CML (near-haploid) | Engineered near-haploid cell line derived from KBM7 (a CML cell line). Unique because most of its genome is haploid — makes genetic screens cleaner since knockouts hit the only copy. |
| BxPC3 | Human | Pancreas | Pancreatic adenocarcinoma | Derived from a 61-year-old female. Adherent. Unlike most pancreatic cancer lines, does NOT carry a KRAS mutation — making it useful for studying KRAS-independent pancreatic cancer biology. |
| HCT116 | Human | Colon | Colorectal carcinoma | Derived from an adult male. Adherent, epithelial. Near-diploid, microsatellite-unstable. Has a KRAS G13D mutation. One of the most commonly used colorectal cancer lines for genetic studies. |
| HEK293T | Human | Kidney (embryonic) | Non-cancerous (transformed) | Human embryonic kidney cells transformed with adenovirus 5 DNA, with an added SV40 large T antigen. Very easy to transfect. The workhorse cell line for protein production, virus packaging, and overexpression experiments. Not truly a kidney cell — likely of neuronal lineage origin. |
| H1 hESC | Human | Blastocyst (embryo) | Non-cancerous (embryonic stem cell) | Human embryonic stem cell line derived from a male blastocyst. Pluripotent — can become any cell type. The cell line used in VCC 2025 challenge. Maintained in undifferentiated state. |
| GBM lines (3 lines) | Human | Brain | Glioblastoma multiforme | Three patient-derived glioblastoma cell lines (specific identities in McFaline-Figueroa paper). GBM is the most aggressive primary brain tumor. |

---

## Dataset-by-Dataset Access Guide

### 1. Replogle et al. 2022 — K562 + RPE1 CRISPRi

**Organism:** Human  
**Cell lines:**
- **K562** — CML blood cancer cells from bone marrow. The default cell line for Perturb-seq work.
- **RPE1** — Normal (non-cancerous) retinal epithelial cells. Provides a non-cancer baseline.

**Format:** h5ad (AnnData), ready for scanpy  
**Total size:** ~11.5 GB across 3 files  
**Ease of access:** Easy — free, no auth, h5ad format, direct download  
**Files on Zenodo (`zenodo.org/records/7041849`):**
- `ReplogleWeissman2022_K562_gwps.h5ad` — **8.8 GB** (genome-wide, ~9,866 gene perturbations in K562)
- `ReplogleWeissman2022_K562_essential.h5ad` — **1.5 GB** (~2,057 essential gene perturbations in K562)
- `ReplogleWeissman2022_rpe1.h5ad` — **1.2 GB** (~2,393 essential gene perturbations in RPE1)

**What's inside:** scRNA-seq count matrices (sparse CSR). Created with scanpy 1.9.1. Preprocessed by scPerturb project. Control cells carry non-targeting guides. Perturbation identity is in obs metadata.

**How to access:**
- Direct download from Zenodo (no auth needed)
- Or via `pertpy.data.replogle_2022_k562_gwps()` (would need to `pip install pertpy`)
- Load with `sc.read_h5ad("ReplogleWeissman2022_K562_gwps.h5ad")`

**Relevance:** The foundational Perturb-seq dataset. K562 and RPE1 are two distinct cell types (cancer blood vs. normal epithelial). The genome-wide K562 screen is the most comprehensive single-cell CRISPRi dataset from this era. Many benchmark methods (GEARS, scGPT) train on this.

---

### 2. Nadig et al. 2025 — HepG2 + Jurkat CRISPRi

**Organism:** Human  
**Cell lines:**
- **HepG2** — Liver cancer cells. Retain liver-specific functions like drug metabolism, making perturbation responses liver-relevant.
- **Jurkat** — T-cell leukemia cells. The standard model for immune/T-cell signaling. Perturbation responses will reflect immune cell biology.

**Format:** Processed single-cell data on GEO (likely h5ad or mtx+barcodes+features)  
**Cell counts:** HepG2: 145,473 cells | Jurkat: 262,956 cells (~408k total)  
**Perturbations:** ~2,393 essential genes per cell line  
**Ease of access:** Medium — free but GEO has CAPTCHA issues for programmatic access, format may need conversion  

**What's inside:** CRISPRi Perturb-seq of essential genes. Aligned with Cell Ranger 4.0.0. Guide assignment via Poisson-Gaussian mixture model. Quality-filtered (removed low UMI and high mito cells). Controls carry non-targeting guide RNAs.

**How to access:**
- GEO accession: **GSE264667**
- Download via GEO's FTP or use `GEOparse` Python package
- Raw data on SRA under BioProject PRJNA1100571
- May need to convert from GEO's format to h5ad if not already provided as such

**Relevance:** Adds HepG2 (liver) and Jurkat (T-cell) — two biologically distinct cell types not in Replogle. Same essential-gene CRISPRi approach. Liver and immune cells have very different gene regulatory programs from blood cancer (K562) or epithelial (RPE1) cells.

---

### 3. Feng et al. 2024 — 34 iPSC lines CRISPRi

**Organism:** Human  
**Cell lines:**
- **34 iPSC lines** — Induced pluripotent stem cells from 26 healthy donors of diverse genetic backgrounds. These are reprogrammed adult cells returned to a pluripotent state. Biologically the closest available model to the H1 hESC line used in VCC 2025 — both are pluripotent human stem cells, just derived differently (iPSC from adult tissue reprogramming vs. hESC from embryos).

**Format:** Count matrices + metadata (NOT h5ad — raw UMI counts, guide counts, cell metadata as separate files)  
**Scale:** 7,226 gene knockdowns across 34 iPSC lines from 26 genetic backgrounds  
**Ease of access:** Medium — free Figshare download, but requires assembling AnnData objects from raw count matrices + metadata  

**What's inside:** Two Perturb-seq screens (genome-wide + targeted) in iPSCs. For each screen, 3 file types: RNA UMI Counts, Guide UMI Counts, Cell Metadata.

**How to access:**
- Count Data: `https://figshare.com/s/bab7c3f17a5fd284fc91`
- Transcriptional Changes: `https://figshare.com/s/14edeeab56eb8a885df3`
- Direct download from Figshare (no auth needed)

**Important note:** Data is in raw count matrix format, not h5ad. Would need to construct AnnData objects from the count matrices and metadata files.

**Relevance:** **Extremely relevant** — iPSCs are the closest cell type to hESCs (VCC 2025 used H1 hESC). 34 lines across 26 donors gives cross-genetic-background generalization data. CRISPRi mechanism matches the VCC challenge. Also uniquely valuable because it shows how the same perturbation plays out differently across genetic backgrounds — exactly the kind of generalization VCC 2026 tests.

---

### 4. Jiang et al. 2025 — 6 cancer cell lines Perturb-seq

**Organism:** Human  
**Cell lines:**
- **A549** — Lung cancer (adenocarcinoma). Epithelial, adherent. Expresses lung markers.
- **MCF7** — Breast cancer (ER+ adenocarcinoma). Estrogen-responsive.
- **HT29** — Colon cancer (adenocarcinoma). Can differentiate into gut-like cells.
- **HAP1** — Near-haploid CML-derived. Unique genetics (one copy of most genes).
- **BxPC3** — Pancreatic cancer (KRAS-wild-type, unusual for pancreatic cancer).
- **K562** — CML blood cancer (same line as Replogle, enabling cross-study comparison).

**Format:** Seurat .rds objects (R format, NOT h5ad)  
**Total size:** ~20.5 GB across 5 pathway-specific files  
**Ease of access:** Medium-hard — free Zenodo download, but requires R-to-Python format conversion  
**Files on Zenodo (`doi.org/10.5281/zenodo.14518762`):**
- IFNB pathway: 4.3 GB
- IFNG pathway: 2.9 GB
- TNFA pathway: 4.7 GB
- TGFB1 pathway: 2.6 GB
- INS pathway: 5.6 GB
- Plus DE results (324 MB zip), bulk RNA-seq (3.4 MB), guide capture protocol (PDF)

**What's inside:** 1,500+ perturbations across all 6 cell lines in 5 signaling pathway contexts (interferon-beta, interferon-gamma, TNF-alpha, TGF-beta1, insulin). Cell metadata includes `guide_identity` (target gene, or "neg" for controls), `cell_type` (cell line), pathway info. Uses their Mixscale framework to weight perturbation efficiency.

**How to access:**
- Download .rds files from Zenodo (no auth)
- **Conversion needed**: Seurat → AnnData. Options:
  - Use `sceasy` R package: `sceasy::convertFormat(obj, from="seurat", to="anndata", outFile="output.h5ad")`
  - Use `SeuratDisk`: write to h5Seurat then convert
  - Or use `rpy2` in Python to load R objects
- Controls: cells with `guide_identity = "neg"`

**Relevance:** The only dataset with 6 different cell lines (spanning lung, breast, colon, blood, pancreas) perturbed with the same set of pathway regulators. Ideal for learning cross-cell-line perturbation responses — the data directly tests whether the same gene knockdown has different effects in different tissues. The R format adds friction but the biological value for cross-cell-line generalization is high.

---

### 5. X-Atlas/Orion (Huang et al. 2025) — HCT116 + HEK293T CRISPRi

**Organism:** Human  
**Cell lines:**
- **HCT116** — Colorectal cancer cells. Near-diploid, microsatellite-unstable, carries KRAS G13D mutation. Widely used in genetic studies because of clean genetics.
- **HEK293T** — Transformed embryonic kidney cells (not truly kidney — likely neuronal lineage). The most commonly used line for transfection and protein production. Extremely well-characterized.

**Format:** h5ad on Figshare, Parquet on HuggingFace  
**Total size:** ~126 GB (HuggingFace Parquet) | h5ad files on Figshare (size TBD, likely 20-40 GB)  
**Scale:** 8M cells, 18,903 gene perturbations (all protein-coding genes), median 16,000 UMIs/cell  
**Ease of access:** Easy — free, no auth, h5ad on Figshare (but very large files)  

**h5ad files on Figshare (`doi.org/10.25452/figshare.plus.29190726`):**
- `HCT116_filtered_dual_guide_cells.h5ad` — cells with 2 sgRNAs targeting the same gene
- `HEK293T_filtered_dual_guide_cells.h5ad` — same for HEK293T
- MD5 checksums provided

**AnnData obs columns:**
- `sample` — GEM batch
- `guide_target` — guide RNA identity
- `gene_target` — gene being targeted
- `n_genes_by_counts`, `total_counts`, `total_counts_mt`, `pct_counts_mt` — QC metrics
- `pass_guide_filter` — boolean, two guides from same pair

**Controls:** 1,026 non-targeting control pairs

**How to access (two options):**
1. **Figshare h5ad** (recommended for our workflow): Direct download, load with `sc.read_h5ad()`
2. **HuggingFace Parquet** (streaming option):
   ```python
   # Would need: pip install datasets huggingface_hub
   from datasets import load_dataset
   ds = load_dataset("Xaira-Therapeutics/X-Atlas-Orion", split="HCT116", streaming=True)
   ```

**Relevance:** The largest public Perturb-seq dataset. Genome-wide CRISPRi in 2 cell lines. Deeply sequenced. HCT116 (colon cancer) and HEK293T (transformed embryonic kidney) are biologically very different from each other and from the other datasets' cell lines. The preferred training resource for VCC 2026. But very large — may need to work with subsets or stream.

---

### 6. Norman et al. 2019 — K562 CRISPRa (Tier 2)

**Organism:** Human  
**Cell lines:**
- **K562** — Same CML blood cancer line as in Replogle. But here perturbed with CRISPRa (gene activation) instead of CRISPRi (gene interference).

**Format:** Available through GEARS/pertpy data loaders, or GEO (GSE133344)  
**Scale:** ~100k cells, 105 single + 131 combinatorial CRISPRa perturbations  
**Ease of access:** Easy if using GEARS/pertpy loaders; otherwise standard GEO download  

**How to access:**
- GEO: GSE133344 (raw/processed)
- `pertpy` data loader (if installed)
- GEARS built-in data loader: `from gears import PertData; pert_data = PertData('./data'); pert_data.load(data_name='norman')`

**Relevance:** Classic benchmark but CRISPRa (activation) ≠ CRISPRi (interference). CRISPRa turns genes UP, CRISPRi turns genes DOWN — different downstream effects. Useful for learning gene-gene interaction structure (which genes co-regulate), less directly applicable to VCC's CRISPRi task.

---

### 7. McFaline-Figueroa et al. 2024 — GBM lines sci-Plex-GxE (Tier 2)

**Organism:** Human  
**Cell lines:**
- **3 patient-derived GBM lines** — Glioblastoma multiforme, the most lethal primary brain tumor. Patient-derived lines better reflect tumor heterogeneity than established lines.

**Format:** Need to check paper data availability (likely GEO)  
**Scale:** 3 GBM cell lines, 522 kinase perturbations x 4 small molecules, ~1M cells  
**Ease of access:** Unknown — need to check paper  

**How to access:** Data availability section of the [Cell Genomics paper](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00339-7). Likely deposited on GEO.

**Relevance:** Combined genetic + chemical perturbation (sci-Plex-GxE). Unique because it measures how genetic background changes drug response — a gene-by-environment interaction. Less directly applicable to VCC's pure CRISPRi task but intellectually interesting for understanding context-dependent perturbation effects.

---

## Practical Considerations

### Storage requirements
| Dataset | Size | Format | Ease |
|---------|------|--------|------|
| Replogle (all 3) | ~11.5 GB | h5ad (ready) | Easy |
| Nadig (HepG2+Jurkat) | ~2-5 GB (est.) | GEO (may need conversion) | Medium |
| Feng (34 iPSC lines) | Unknown (est. 5-15 GB) | Count matrices (need to build AnnData) | Medium |
| Jiang (6 cell lines) | ~20.5 GB | Seurat .rds (need R→Python conversion) | Medium-hard |
| X-Atlas (2 cell lines) | ~20-40 GB (h5ad) | h5ad (ready) | Easy (but large) |
| **Total estimate** | **~60-90 GB** | | |

### Format compatibility
- **Ready to use** (h5ad): Replogle, X-Atlas/Orion
- **Needs conversion from R**: Jiang (Seurat .rds → h5ad)
- **Needs assembly**: Feng (count matrices + metadata → AnnData)
- **Needs format check**: Nadig (GEO download, likely mtx or h5ad)

### New packages needed
To access all datasets, would need to add:
- `huggingface_hub` / `datasets` — for X-Atlas streaming (optional if using Figshare h5ad)
- `pertpy` — for convenient Replogle/Norman loading (optional if downloading from Zenodo)
- R + `sceasy` or `SeuratDisk` — for converting Jiang Seurat objects (one-time conversion)
