r"""
Master Script for Generating All Figures (Figures 1, 2, 3, and 4).

Generates all publication-ready PDF figures for the manuscript:
- Figure 1: The smallest conflicting pair of trees (N = 7).
- Figure 2: (W(T), a(T)) scatter plot for all 106 trees of order 10.
- Figure 3: Algebraic connectivity curves for Theorems 1 and 2 families.
- Figure 4: Complete diameter-by-diameter classification grid (N \le 19).

All figures are saved to the 'figures/' directory.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig1 import generate_figure_1
from fig2 import generate_figure_2
from fig3 import generate_figure_3
from fig4 import generate_figure_4
from utils import FIGURES_DIR, spider  # re-export spider for backwards compatibility


def generate_all_figures():
    """Run all individual figure generation routines."""
    print("=" * 65)
    print("GENERATING MANUSCRIPT FIGURES INTO 'figures/'")
    print("=" * 65)
    
    print("\n[1/4] Generating Figure 1 (Smallest conflicting pair, N = 7)...")
    generate_figure_1(FIGURES_DIR / "fig1_smallest_conflict.pdf")
    
    print("\n[2/4] Generating Figure 2 (Scatter plot for order N = 10)...")
    generate_figure_2(FIGURES_DIR / "fig2_scatter_n10.pdf")
    
    print("\n[3/4] Generating Figure 3 (Theorems 1 and 2 spectral curves)...")
    generate_figure_3(FIGURES_DIR / "fig3_family_an_bn.pdf")
    
    print("\n[4/4] Generating Figure 4 (Diameter-class classification grid)...")
    generate_figure_4(output_path=FIGURES_DIR / "fig4_diameter_grid.pdf")
    
    print("\n" + "=" * 65)
    print("ALL FIGURES SUCCESSFULLY GENERATED!")
    print("=" * 65)


if __name__ == "__main__":
    generate_all_figures()
