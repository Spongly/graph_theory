r"""
Figure 2: Scatter Plot of (W(T), a(T)) for Trees of Order N = 10.

Plots algebraic connectivity a(T) versus Wiener index W(T) for all 106
non-isomorphic trees of order 10, coloured by diameter d \in {2, ..., 9}.
Output file: figures/fig2_scatter_n10.pdf
"""

from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import FIGURES_DIR, compute_metrics


def generate_figure_2(output_path=None):
    """Generate and save Figure 2."""
    if output_path is None:
        output_path = FIGURES_DIR / "fig2_scatter_n10.pdf"
        
    plt.rcParams.update({
        "font.size": 9,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    
    # Compute metrics for all 106 trees of order 10
    tree_records = []
    for tree in nx.nonisomorphic_trees(10):
        a_val, w_val = compute_metrics(tree)
        diam = nx.diameter(tree)
        tree_records.append((a_val, w_val, diam))
        
    fig, ax = plt.subplots(figsize=(3.3, 2.8))
    cmap = plt.get_cmap("viridis")
    
    # Plot points grouped by diameter d
    for d in range(2, 10):
        group = [rec for rec in tree_records if rec[2] == d]
        if not group:
            continue
        w_coords = [rec[1] for rec in group]
        a_coords = [rec[0] for rec in group]
        
        ax.scatter(
            w_coords,
            a_coords,
            s=13,
            color=cmap((d - 2) / 7.0),
            label=f"$d={d}$",
            zorder=3,
            edgecolors="none",
        )
        
    ax.set_xlabel(r"Wiener index $W(T)$")
    ax.set_ylabel(r"algebraic connectivity $a(T)$")
    ax.legend(frameon=False, fontsize=6.5, ncol=2, loc="upper right")
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved Figure 2 to {output_path}")


if __name__ == "__main__":
    generate_figure_2()
