r"""
Figure 4: Classification of Diameter Classes for Orders 7 \le N \le 19.

Visualizes the diameter-by-diameter classification of tree robustness orderings:
- Blue circles: Agreement between a-ordering and W-ordering (proved for d \in {2, 3, N-1}).
- Dark orange squares: Conflict proved theoretically (Theorems 1 and 2).
- Light orange squares: Conflict computer-verified only.

Reads data from: code/scan.json
Output file: figures/fig4_diameter_grid.pdf
"""

import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import CODE_DIR, FIGURES_DIR


def is_theoretically_proved(N: int, d: int) -> bool:
    """
    Check if the conflict at (N, d) is covered by Theorem 1 (even d) or Theorem 2 (odd d).
    """
    if d % 2 == 0:
        return N >= 3 * d // 2 + 1
    else:
        return N >= 2 * d


def generate_figure_4(json_path=None, output_path=None):
    """Generate and save Figure 4 from scan.json data."""
    if json_path is None:
        json_path = CODE_DIR / "scan.json"
    if output_path is None:
        output_path = FIGURES_DIR / "fig4_diameter_grid.pdf"
        
    if not json_path.exists():
        raise FileNotFoundError(f"Missing required data file: {json_path}")
        
    with open(json_path, "r", encoding="utf-8") as f:
        scan_data = json.load(f)
        
    plt.rcParams.update({
        "font.size": 8.5,
        "font.family": "serif",
        "mathtext.fontset": "cm",
    })
    
    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    
    categories = {
        "agreement": {
            "color": "#4a6fa5",
            "marker": "o",
            "label": r"agreement (computed; proved for $d\in\{2,3,N-1\}$)",
        },
        "proved_conflict": {
            "color": "#d1692a",
            "marker": "s",
            "label": r"conflict, proved (Theorems 1 and 2)",
        },
        "unproved_conflict": {
            "color": "#f2b98f",
            "marker": "s",
            "label": r"conflict, computer-verified only",
        },
    }
    
    used_labels = set()
    
    for str_N, row in scan_data.items():
        N = int(str_N)
        for str_d, record in row.items():
            d = int(str_d)
            
            if not record["conflict"]:
                cat_key = "agreement"
            else:
                cat_key = "proved_conflict" if is_theoretically_proved(N, d) else "unproved_conflict"
                
            cfg = categories[cat_key]
            lbl = cfg["label"] if cat_key not in used_labels else None
            used_labels.add(cat_key)
            
            edge_color = "k" if cat_key != "unproved_conflict" else "#999"
            ax.scatter(
                N,
                d,
                s=26,
                c=cfg["color"],
                marker=cfg["marker"],
                edgecolors=edge_color,
                linewidths=0.4,
                label=lbl,
            )
            
    ax.set_xlabel(r"order $N$")
    ax.set_ylabel(r"diameter $d$")
    ax.set_xticks(range(7, 20))
    ax.set_yticks(range(2, 19, 2))
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.grid(alpha=0.2, lw=0.4)
    
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
        
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved Figure 4 to {output_path}")


if __name__ == "__main__":
    generate_figure_4()
