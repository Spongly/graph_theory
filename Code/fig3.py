r"""
Figure 3: Algebraic Connectivity of Trees A and B in Theorems 1 and 2.

Plots a(A) and a(B) as functions of order N for:
- Left panel (Theorem 1, even diameter d \in {4, 6, 8}):
    A = T(h, h, h) meets the plateau \lambda_h (solid line).
    B = T(h+1, h-1) strictly dips below the plateau (dashed line).

- Right panel (Theorem 2, odd diameter d \in {5, 7, 9}):
    A = T(h+1, h, h, h) satisfies a(A) >= \lambda_{h+1} (solid line).
    B = T(h+2, h-1, h-1) strictly satisfies a(B) < \lambda_{h+1} (dashed line).
    Dotted horizontal line denotes the threshold \lambda_{h+1}.

Output file: figures/fig3_family_an_bn.pdf
"""

from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import FIGURES_DIR, spider, compute_algebraic_connectivity, pinned_path_eigenvalue


def generate_figure_3(output_path=None):
    """Generate and save Figure 3."""
    if output_path is None:
        output_path = FIGURES_DIR / "fig3_family_an_bn.pdf"
        
    plt.rcParams.update({
        "font.size": 9.5,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    
    color_map = {2: "#4a6fa5", 3: "#d1692a", 4: "#3a8f5c"}
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8), sharey=False)
    
    # Left panel: Theorem 1 (Even diameter d = 2h)
    for h in (2, 3, 4):
        orders = list(range(3 * h + 1, 41))
        a_A = [compute_algebraic_connectivity(spider([h, h, h], N)) for N in orders]
        a_B = [compute_algebraic_connectivity(spider([h + 1, h - 1], N)) for N in orders]
        
        axes[0].plot(orders, a_A, "-", c=color_map[h], lw=1.8, label=f"$d={2*h}$")
        axes[0].plot(orders, a_B, "--", c=color_map[h], lw=1.4)
        
    # Right panel: Theorem 2 (Odd diameter d = 2h + 1)
    for h in (2, 3, 4):
        orders = list(range(4 * h + 2, 46))
        a_A = [compute_algebraic_connectivity(spider([h + 1, h, h, h], N)) for N in orders]
        a_B = [compute_algebraic_connectivity(spider([h + 2, h - 1, h - 1], N)) for N in orders]
        threshold = pinned_path_eigenvalue(h + 1)
        
        axes[1].plot(orders, a_A, "-", c=color_map[h], lw=1.8, label=f"$d={2*h+1}$")
        axes[1].plot(orders, a_B, "--", c=color_map[h], lw=1.4)
        axes[1].axhline(threshold, c=color_map[h], lw=0.8, ls=":")
        
    axes[0].set_title(r"Theorem 1: $A$ solid, $B$ dashed", fontsize=9)
    axes[1].set_title(r"Theorem 2 (dotted: $\lambda_{h+1}$)", fontsize=9)
    
    axes[0].set_ylim(0.07, 0.44)
    axes[1].set_ylim(0.06, 0.29)
    
    for ax in axes:
        ax.set_xlabel("order $N$")
        ax.legend(frameon=False, fontsize=8.5, ncol=3, loc="upper right", columnspacing=1)
        
    axes[0].set_ylabel("algebraic connectivity")
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved Figure 3 to {output_path}")


if __name__ == "__main__":
    generate_figure_3()
