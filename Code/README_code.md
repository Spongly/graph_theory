# Computational Pipeline & Scripts

This folder contains Python scripts accompanying the paper:
**"Algebraic Connectivity Versus the Kirchhoff Index: When Do They Rank Trees Differently?"**

## Prerequisites

- **Python**: 3.10+
- **Packages**:
  ```bash
  pip install numpy scipy sympy networkx matplotlib
  ```

---

## File Overview & Mapping to Manuscript

| Script | Paper Section / Result | Description |
|---|---|---|
| [`utils.py`](utils.py) | Shared | Common helpers for tree construction (`spider`, `comet`), spectral computations, Wiener indices, and path resolution. |
| [`verify.py`](verify.py) | **Theorems 1 & 2** | Verifies the explicit conflicting families of trees across $h = 2..9$ and orders $N$, and checks the closed-form secular function identity $\Phi_B(\lambda_h) = -\ell \frac{\lambda_h}{1 - \lambda_h}$. |
| [`verify_certificates.py`](verify_certificates.py) | **Propositions 2 & 4** | Certifies conflicts in exact rational arithmetic via the inertia lemma for `witnesses.csv` (all 91 classes $N \le 19$) and `window_certificates.csv` (1,859 pairs). |
| [`generate_witnesses.py`](generate_witnesses.py) | **Proposition 4 / Table** | Searches and produces `witnesses.csv` with exact rational separators for all $N \le 19, 4 \le d \le N-3$. |
| [`generate_window_certificates.py`](generate_window_certificates.py) | **Proposition 2** | Constructs and certifies `window_certificates.csv` for the odd-diameter window ($h=3..60$, $r=2..h$ and $h=61..150$, $r=2$). |
| [`gaps.py`](gaps.py) | **Theorems 1 & 2** | Symbolically derives Wiener index gaps $\Delta W_A$ and $\Delta W_B$ using SymPy and asserts exact equality with NetworkX calculations. |
| [`tables.py`](tables.py) | **Tables 1 & 2** | Exhaustively enumerates small connected graphs ($n \le 7$) and non-isomorphic trees ($7 \le N \le 14$), computing conflict percentages. |
| [`scan.py`](scan.py) | **Section 4 & Figure 4** | Scans trees by diameter in $O(n \log n)$ time to test whether diameter classes contain conflicts, saving records to `scan.json`. |
| [`cover.py`](cover.py) | **Corollary 1** | Compares empirical diameter conflicts with theoretical bounds: $N \ge 3d/2 + 1$ (even $d$) and $N \ge 2d$ (odd $d$). |
| [`comet.py`](comet.py) | **Open Problem 1 (d = N - 2)** | Tests comet trees $P_{N-1}$ with a pendant leaf for $N = 7..120$, verifying that $a$ strictly increases while $W$ strictly decreases toward the center. |
| [`near.py`](near.py) | **Open Problem 2 (d = N - 3)** | Analyzes large-diameter trees with $d = N - 3$ for $N \in \{14, 16, 18, 20\}$ and identifies pairs with minimal $\Delta W$. |
| [`fig1.py`](fig1.py) | **Figure 1** | Generates `figures/fig1_smallest_conflict.pdf` (smallest conflicting pair, $N = 7$). |
| [`fig2.py`](fig2.py) | **Figure 2** | Generates `figures/fig2_scatter_n10.pdf` (scatter plot of $(W(T), a(T))$ for all 106 trees of order 10). |
| [`fig3.py`](fig3.py) | **Figure 3** | Generates `figures/fig3_family_an_bn.pdf` (algebraic connectivity vs. $N$ for Theorems 1 and 2). |
| [`fig4.py`](fig4.py) | **Figure 4** | Generates `figures/fig4_diameter_grid.pdf` (diameter classification grid $N \le 19$). |
| [`figs.py`](figs.py) | **Figures 1–4** | Master generator that builds all figures sequentially. |

---

## Running the Code

All scripts can be executed directly from the project root or from inside `code/`:

```bash
# 1. Verify Theorems 1 & 2 and secular equations
python code/verify.py

# 2. Check symbolic Wiener index gap formulas
python code/gaps.py

# 3. Reproduce Tables 1 and 2
python code/tables.py

# 4. Generate all figures
python code/figs.py

# 5. Check comet trees (d = N - 2)
python code/comet.py
```
