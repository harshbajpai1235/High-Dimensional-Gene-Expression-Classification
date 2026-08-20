# High-Dimensional Gene Expression Classification

MATH 748 project on **regularization, dimension reduction, and model comparison** using the NCI60 cancer cell-line gene expression data. The setting is classic high-dimension, small-sample learning:

$$p \gg n \quad (6{,}830 \text{ genes},\ 64 \text{ samples})$$

The work compares unregularized logistic regression, Ridge, Lasso, and PCA + logistic regression for cancer-type classification, then extends the best approach from a binary task to an 8-class problem.

**Author:** Harsh Bajpai

---

## Question

Can cancer type be predicted from gene expression when there are far more genes than samples, and which regularization strategy actually helps?

The analysis studies:

1. **Supervised classification** of cancer type from expression profiles
2. **High-dimensional modeling** — Ridge, Lasso, Elastic-style penalties, PCA, and the bias–variance tradeoff in $\lambda$ / $C$
3. **Stability and interpretability** — sparse gene selection vs shrinking all coefficients vs compressing all genes

---

## Dataset

[NCI60](https://islp.readthedocs.io/) via `ISLP.load_data("NCI60")`:

| | |
| --- | --- |
| Samples $n$ | 64 cancer cell lines |
| Features $p$ | 6,830 gene expression measurements |
| Labels | 14 cancer types (highly imbalanced) |

Rare classes (1–2 samples) make full 14-way classification unstable. The notebook therefore:

- Starts with a **balanced binary** task: **RENAL vs NSCLC** (9 + 9 = 18 samples)
- Then a **filtered multi-class** task: 8 types with at least 5 samples each (**57 samples**)

No near-constant genes were found (minimum gene variance ≈ 0.04). Features are standardized to mean 0 and variance 1 before regularized models and PCA.

PCA on the full matrix: PC1 ≈ **11.4%** of variance; first 10 PCs ≈ **46%**. MELANOMA and LEUKEMIA separate somewhat in the PC1–PC2 plane; RENAL, NSCLC, BREAST, and COLON overlap.

---

## Methods

| Stage | Approach |
| --- | --- |
| Baseline | Unregularized logistic regression, 5-fold CV |
| L2 | Ridge logistic regression, grid over $C$ |
| L1 | Lasso logistic regression (`liblinear`), grid over $C$, sparse gene counts |
| Dimension reduction | PCA (2 / 5 / 10 components) + logistic regression |
| Multi-class | Multinomial / one-vs-rest logistic and Lasso on 8 classes |
| Diagnostics | Confusion matrix, per-class precision / recall / F1 |

---

## Key results

### Binary: RENAL vs NSCLC (18 × 6,830)

| Model | Best CV accuracy |
| --- | --- |
| Logistic (no penalty) | ~73% (unstable folds, 50–100%) |
| Ridge (tuned) | ~73% |
| PCA + logistic | ~73% (2–10 components) |
| **Lasso (tuned, $C = 10$)** | **~85%** |

Lasso at $C = 10$ keeps **54 genes** (more than 99% of features set to zero). Stronger penalties underfit (~50%); weaker penalties select more genes and lose a bit of accuracy. Ridge shrinks coefficients but keeps all 6,830 genes, so noise stays in the model. PCA compresses rather than drops irrelevant genes, so it matches Ridge rather than Lasso.

Selected gene *sets* can shift slightly across refits — expected when $n = 18$ and genes are correlated. Treat the sparse list as illustrative, not a unique biomarker panel.

### Multi-class: 8 cancer types (57 × 6,830)

| Model | CV accuracy |
| --- | --- |
| Multinomial logistic | ~67% (folds ~45–91%) |
| **Lasso (best $C = 100$)** | **~65%** |

Multi-class needs weaker regularization (larger $C$) because more genes are required to separate eight boundaries. Performance is uneven: COLON and LEUKEMIA are nearly separable; MELANOMA, OVARIAN, and NSCLC are moderate; BREAST and CNS are hard — consistent with overlap in expression space.

**Takeaway:** in this $p \gg n$ setting, **feature selection (Lasso) beats coefficient shrinkage (Ridge) and compression (PCA)**. Binary tasks look strong; realistic multi-class cancer typing is much harder with $n \approx 60$.

---

## Repository contents

| File | Description |
| --- | --- |
| `Math_748_Project_Harsh_Bajpai.ipynb` | Full analysis notebook |
| `Math_748_HB_Final_Project_Report.pdf` | Written report |
| `Math_748_Project_Harsh_Bajpai.pdf` | Exported notebook PDF |
| `Harsh Bajpai Math 748 Project Presentation.pptx` | Project slides |
| `Project_Proposal_HB.pdf` / `.docx` | Original project proposal |

The NCI60 matrix is **not** stored in the repo; the notebook downloads it through ISLP.

---

## How to run

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Open `Math_748_Project_Harsh_Bajpai.ipynb` in Jupyter, VS Code, or Cursor. The first data cell is:

```python
from ISLP import load_data
nci = load_data("NCI60")
```

---

## Stack

Python, NumPy, pandas, matplotlib, scikit-learn, ISLP (NCI60).

---

## Disclaimer

Academic analysis of a public cell-line expression dataset. It is not a clinical diagnostic model and selected genes should not be treated as validated biomarkers.
