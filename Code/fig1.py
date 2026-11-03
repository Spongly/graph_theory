r"""
Figure 1: The Smallest Conflicting Pair of Trees (Order N = 7).

Draws the two conflicting trees of order 7:
- Left panel: Spider tree A_7 = T(2, 2, 2)
    Algebraic connectivity a \approx 0.3820, Wiener index W = 48.
- Right panel: Tree B_7 = T(3, 1, 1, 1)
    Algebraic connectivity a \approx 0.2955, Wiener index W = 46.

Here a(A_7) > a(B_7) while W(A_7) > W(B_7), demonstrating a conflict.
Output file: figures/fig1_smallest_conflict.pdf
"""

from math import cos, pi, sin
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import FIGURES_DIR, spider, compute_metrics


def draw_tree_panel(ax, pos, G, title):
    """Draw network nodes and edges on the specified matplotlib axis."""
    nx.draw_networkx_edges(G, pos, ax=ax, width=1.4)
    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        node_size=110,
        node_color="#f2f2f2",
        edgecolors="k",
        linewidths=1.1,
    )
    ax.set_title(title, fontsize=11)
    ax.axis("off")
    ax.set_aspect("equal")


def generate_figure_1(output_path=None):
    """Generate and save Figure 1."""
    if output_path is None:
        output_path = FIGURES_DIR / "fig1_smallest_conflict.pdf"
        
    plt.rcParams.update({
        "font.size": 9,
        "font.family": "serif",
        "mathtext.fontset": "cm",
    })
    
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.6))
    
    # Tree A_7: T(2, 2, 2)
    tree_A = spider([2, 2, 2], 7)
    pos_A = {0: (0.0, 0.0)}
    # 3 symmetric legs at angles 90, 210, 330 degrees
    for leg_idx, angle_deg in enumerate([90, 210, 330]):
        angle_rad = angle_deg * pi / 180.0
        for step in (1, 2):
            node_id = 1 + 2 * leg_idx + (step - 1)
            pos_A[node_id] = (step * cos(angle_rad), step * sin(angle_rad))
            
    a_A, w_A = compute_metrics(tree_A)
    draw_tree_panel(
        axes[0],
        pos_A,
        tree_A,
        rf"$A_7 = T(2,2,2)$: $a \approx {a_A:.4f},\ W = {w_A}$",
    )
    
    # Tree B_7: T(3, 1, 1, 1)
    tree_B = spider([3, 1, 1, 1], 7)
    pos_B = {
        0: (0.0, 0.0),
        1: (1.0, 0.0),
        2: (2.0, 0.0),
        3: (3.0, 0.0),
        4: (-0.8, 0.9),
        5: (-1.1, 0.0),
        6: (-0.8, -0.9),
    }
    a_B, w_B = compute_metrics(tree_B)
    draw_tree_panel(
        axes[1],
        pos_B,
        tree_B,
        rf"$B_7 = T(3,1,1,1)$: $a \approx {a_B:.4f},\ W = {w_B}$",
    )
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved Figure 1 to {output_path}")


if __name__ == "__main__":
    generate_figure_1()
