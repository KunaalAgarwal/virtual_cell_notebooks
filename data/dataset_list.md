# List of datasets

As we get more datasets add to the list so we know what we are using and potentially paying for.

## Virtual cell challenge validation set 2026 - Arc Institute
- Access (already in the repo but for principles sake): ```vcc datasets list```, ```vcc datasets download controls```

## Virtual cell atlas - Arc Institute
- Be careful about pulling the entire thing as it will exceed what we are able to freely access at once.
- https://console.cloud.google.com/marketplace/product/bigquery-public-data/arc-institute?project=gcp-public-data-arc-institute&pli=1
- https://github.com/ArcInstitute/arc-virtual-cell-atlas
- Includes the 2025 vcc ([documentation](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/virtual-cell-challenge/README.md)), scBaseCount ([documentation](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/scBaseCount/README.md)), and Tahoe-100M datasets ([documentation](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/tahoe-100M/README.md))

## Public perturbation datasets

Listed on the VCC datasets page: https://virtualcellchallenge.org/datasets

Papers also stored in `virtual_cell_background/data/` 

### CRISPRi, multi-cell-line, large-scale

**Replogle et al. 2022** — "Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq" (Cell)
- Cell lines: K562 (genome-wide ~9,866 genes), K562 (essential ~2,057 genes), RPE1 (essential ~2,393 genes)
- Perturbation type: CRISPRi
- Data: Zenodo `zenodo.org/records/7041849`, GEO GSE146194
- Also available via `pertpy.data.replogle_2022_k562_gwps()`

**Nadig et al. 2025** — "Transcriptome-wide characterization of genetic perturbations" (Nature Genetics)
- Cell lines: HepG2 (~2,393 essential gene perturbations), Jurkat (~2,393 essential gene perturbations)
- Perturbation type: CRISPRi
- Data: GEO GSE264667

**Feng et al. 2024** — "A genome-scale single cell CRISPRi map of trans gene regulation across human pluripotent stem cell lines" (Cell Genomics 2025)
- Cell lines: 34 iPSC lines from 26 genetic backgrounds
- Perturbation type: CRISPRi, 7,226 genes, 20,000+ guide RNAs
- Data: Figshare `figshare.com/s/bab7c3f17a5fd284fc91` (counts), `figshare.com/s/14edeeab56eb8a885df3` (transcriptional changes)
- bioRxiv: https://www.biorxiv.org/content/10.1101/2024.11.28.625833v1

**Jiang et al. 2025** — "Systematic reconstruction of molecular pathway signatures using scalable single-cell perturbation screens" (Nature Cell Biology)
- Cell lines: A549, MCF7, HT29, HAP1, BxPC3, K562 (6 cell lines)
- Perturbation type: Perturb-seq, 1,500+ perturbations, 5 signaling contexts
- Data: Zenodo `doi.org/10.5281/zenodo.14518762`

**X-Atlas/Orion — Huang et al. 2025** — "Genome-wide Perturb-seq Datasets via a Scalable Fix-Cryopreserve Platform" (bioRxiv)
- Cell lines: HCT116, HEK293T
- Perturbation type: CRISPRi, 18,903 gene perturbations, 8M cells, median 140+ cells/perturbation
- Data: HuggingFace `Xaira-Therapeutics/X-Atlas-Orion`, Figshare
- Largest publicly available Perturb-seq dataset to date
- bioRxiv: https://www.biorxiv.org/content/10.1101/2025.06.11.659105v1

### Useful supplementary

**Norman et al. 2019** — "Exploring genetic interaction manifolds constructed from rich single-cell phenotypes" (Science)
- Cell lines: K562
- Perturbation type: CRISPRa (activation, not interference), 105 single + 131 combinatorial perturbations
- Data: GEO GSE133344
- Useful for gene-gene interaction learning (CRISPRa != CRISPRi)

**McFaline-Figueroa et al. 2024** — "Multiplex single-cell chemical genomics reveals the kinase dependence of the response to targeted therapy" (Cell Genomics)
- Cell lines: 3 glioblastoma cell lines
- Perturbation type: sci-Plex-GxE (combined genetic + chemical), 522 kinase perturbations x 4 small molecules, ~1M cells
- Data: check paper data availability section

### Databases

**PerturBase** — 122 datasets from 46 studies, 24,254 genetic + 230 chemical perturbations, ~5M cells
- Website: http://www.perturbase.cn/
- Paper: https://academic.oup.com/nar/article/53/D1/D1099/7815638

### Summary of unique cell lines across all datasets

| Cell Line | Source Paper | Tissue Origin |
|-----------|-------------|---------------|
| K562 | Replogle 2022, Jiang 2025, Norman 2019 | Chronic myelogenous leukemia |
| RPE1 | Replogle 2022 | Retinal pigmented epithelium |
| HepG2 | Nadig 2025 | Hepatocellular carcinoma |
| Jurkat | Nadig 2025 | T-cell leukemia |
| 34 iPSC lines | Feng 2024 | Induced pluripotent stem cells |
| A549 | Jiang 2025 | Lung adenocarcinoma |
| MCF7 | Jiang 2025 | Breast adenocarcinoma |
| HT29 | Jiang 2025 | Colorectal adenocarcinoma |
| HAP1 | Jiang 2025 | Chronic myelogenous leukemia (near-haploid) |
| BxPC3 | Jiang 2025 | Pancreatic adenocarcinoma |
| HCT116 | X-Atlas 2025 | Colorectal carcinoma |
| HEK293T | X-Atlas 2025 | Embryonic kidney |
| H1 hESC | VCC 2025 | Human embryonic stem cells |
| 3 GBM lines | McFaline-Figueroa 2024 | Glioblastoma |
