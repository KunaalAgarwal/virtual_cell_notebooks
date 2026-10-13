# public_screens/

Borrow known knockdown effects, weighted by how alike two genes look in public Perturb-seq screens. Tested on the 2025 H1 data; results in [RESULTS.md](RESULTS.md).

These are evaluation scripts on H1, not yet a `predict()` model for the 2026 contexts.

| File | What it does |
|---|---|
| `screen_transfer_h1.py` | Leave-one-out over the 150 H1 knockdowns: predict each from the others, weighted by screen similarity, and score with pseudobulk Pearson. |
| `celleval_screen_h1.py` | Builds cell-eval submissions (control cells shifted by the predicted change, with the knockdown-specific part optionally amplified) and scores them on 40 real H1 knockdowns next to the baselines, a shuffled null and the oracle. |
| `screen_bulk_from_cells.py` | Turns a single-cell screen (e.g. the support set's `hepg2.h5`, `jurkat.h5`) into z-scored pseudobulk profiles in Replogle's layout. |
| `h1_pseudobulk.py` | Batch-matched H1 pseudobulk effects (150 knockdowns). Handles the support set's log1p-count file. |
| `h1_extract_perts.py` | Pulls the 40-knockdown H1 test set for cell-eval. |
| `h1_eval_utils.py` | Shared scoring and submission-building helpers. |

The whole pipeline, including downloads, runs in Colab via [`notebooks/public_screens_colab.ipynb`](../../notebooks/public_screens_colab.ipynb).

## Data (all public, no VCC key)

- Replogle 2022 K562 genome-wide and RPE1 pseudobulk: figshare files [35773217](https://ndownloader.figshare.com/files/35773217) and [35775512](https://ndownloader.figshare.com/files/35775512).
- 2025 VCC support set (H1 training cells, Nadig HepG2/Jurkat subsets): `https://storage.googleapis.com/vcc_data_prod/datasets/state/competition_support_set.zip` (8.7 GB).

## Run

```bash
python models/public_screens/h1_pseudobulk.py competition_train.h5 h1_pseudobulk.npz
python models/public_screens/h1_extract_perts.py competition_train.h5 h1_pseudobulk.npz h1_real_test.h5ad
python models/public_screens/screen_bulk_from_cells.py hepg2.h5 hepg2_bulk.h5ad
python models/public_screens/screen_bulk_from_cells.py jurkat.h5 jurkat_bulk.h5ad

python models/public_screens/screen_transfer_h1.py --pseudobulk h1_pseudobulk.npz --genes-from score_genes.csv \
    --screen k562=K562_gwps_normalized_bulk_01.h5ad rpe1=rpe1_normalized_bulk_01.h5ad \
             hepg2=hepg2_bulk.h5ad jurkat=jurkat_bulk.h5ad --out results/

python models/public_screens/celleval_screen_h1.py --real h1_real_test.h5ad --pseudobulk h1_pseudobulk.npz \
    --screens K562_gwps_normalized_bulk_01.h5ad rpe1_normalized_bulk_01.h5ad hepg2_bulk.h5ad jurkat_bulk.h5ad \
    --amp 4 8 --shuffled-amp 8 --out celleval/
```

`score_genes.csv` (column `target`) lists the genes to score on. The reported numbers use the 10,993 genes of the H1 control-cell network; the Colab notebook writes that list.
