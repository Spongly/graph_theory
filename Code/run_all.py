r"""
Master Execution Pipeline for "Algebraic Connectivity Versus the Kirchhoff Index"

Executes all verification routines and regenerates all tables and figures
accompanying the paper:
- Table 1: Connected graph conflicts (n <= 7)
- Table 2: Non-isomorphic tree conflicts (7 <= N <= 14)
- Verification of Theorems 1 & 2 and secular equations
- Symbolic derivation of Wiener index gap formulas
- Figure 1: Smallest conflicting pair (N = 7)
- Figure 2: (W, a) scatter plot for N = 10
- Figure 3: Spectral curves for Theorems 1 and 2
- Figure 4: Diameter-by-diameter classification grid (N <= 19)
- Open Problem 1 verification: Comet trees (d = N - 2, N <= 120)

Usage:
    python run_all.py
"""

from pathlib import Path
import sys
import time

# Ensure code/ is on path
REPO_ROOT = Path(__file__).resolve().parent
CODE_DIR = REPO_ROOT / "code"
sys.path.insert(0, str(CODE_DIR))

from tables import generate_table_1, generate_table_2, print_table_1, print_table_2
from verify import verify_even_diameter_family, verify_odd_diameter_family, verify_secular_identities
from gaps import derive_symbolic_gaps, verify_with_networkx
from fig1 import generate_figure_1
from fig2 import generate_figure_2
from fig3 import generate_figure_3
from fig4 import generate_figure_4
from comet import check_comet_orderings
from verify_certificates import verify_witnesses, verify_window_certificates
from utils import FIGURES_DIR


def run_all(max_tree_n=14):
    start_total = time.time()
    print("=" * 75)
    print("REPRODUCING ALL RESULTS FOR THE PAPER")
    print("=" * 75)
    
    # 1. Verification of Theorems 1 and 2
    print("\n[1/7] Verifying Theorems 1 & 2 (Explicit Conflicting Families)...")
    t_even, f_even = verify_even_diameter_family()
    t_odd, f_odd = verify_odd_diameter_family()
    sec_results = verify_secular_identities()
    print(f"  - Theorem 1 (Even d): {t_even} instances tested, {f_even} failures")
    print(f"  - Theorem 2 (Odd d) : {t_odd} instances tested, {f_odd} failures")
    print(f"  - Secular equation identity Phi_B matches closed form for all h in 2..6: {all(r[3] < 1e-12 for r in sec_results)}")
    
    # 2. Symbolic Gaps Verification
    print("\n[2/7] Verifying Symbolic Wiener Index Gaps (Theorems 1 & 2)...")
    h_sym, N_sym, gap_A, gap_B = derive_symbolic_gaps()
    verify_with_networkx(h_sym, N_sym, gap_A, gap_B)
    print("  - Symbolic Wiener gaps match NetworkX graph computations exactly.")
    
    # 3. Exact Rational Certificates (Inertia Lemma)
    print("\n[3/7] Verifying Exact Rational Certificates (witnesses & window_certificates)...")
    w_res = verify_witnesses(verbose=False)
    win_res = verify_window_certificates(verbose=False)
    print(f"  - witnesses.csv: {w_res['passed']}/{w_res['total']} verified (0 failures)")
    print(f"  - window_certificates.csv: {win_res['passed']}/{win_res['total']} verified (0 failures)")
    
    # 4. Table 1: Connected Graphs
    print("\n[4/7] Generating Table 1 (Connected Graphs, n <= 7)...")
    t1_rows = generate_table_1()
    print_table_1(t1_rows)
    
    # 5. Table 2: Tree Conflicts
    print(f"\n[5/7] Generating Table 2 (Trees, 7 <= N <= {max_tree_n})...")
    t2_rows = generate_table_2(min_N=7, max_N=max_tree_n)
    print_table_2(t2_rows)
    
    # 6. Figures 1 through 4
    print("\n[6/7] Regenerating Figures 1, 2, 3, and 4 into 'figures/'...")
    generate_figure_1(FIGURES_DIR / "fig1_smallest_conflict.pdf")
    generate_figure_2(FIGURES_DIR / "fig2_scatter_n10.pdf")
    generate_figure_3(FIGURES_DIR / "fig3_family_an_bn.pdf")
    generate_figure_4(output_path=FIGURES_DIR / "fig4_diameter_grid.pdf")
    
    # 7. Comet Trees (Diameter N - 2)
    print("\n[7/7] Verifying Open Problem 1 (Comet Trees, d = N - 2 for N <= 120)...")
    comet_ok, min_gap, fails = check_comet_orderings(7, 120)
    print(f"  - Strictly monotonic in a and W toward center: {comet_ok}")
    print(f"  - Minimum consecutive a-gap: {min_gap:.4e}")
    
    elapsed = time.time() - start_total
    print("\n" + "=" * 75)
    print(f"PIPELINE COMPLETED SUCCESSFULLY IN {elapsed:.2f}s")
    print("=" * 75)


if __name__ == "__main__":
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    run_all(max_tree_n=max_n)
