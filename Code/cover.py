r"""
Comparison of Empirical Diameter-Class Conflicts vs. Theoretical Coverage.

For orders N = 7 to 16 and intermediate diameters 4 <= d <= N - 3:
- Empirically detects if any conflicting pair exists within the diameter class.
- Checks whether the diameter class falls within the theoretical coverage:
    - Even d: N >= 3d/2 + 1 (Theorem 1)
    - Odd d:  N >= 2d       (Theorem 2)

This illustrates the 'dark orange' vs 'light orange' regions shown in Figure 4.
"""

from collections import defaultdict
from pathlib import Path
import sys
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import compute_metrics, has_conflict


def is_covered_by_theorems(N: int, d: int) -> bool:
    """
    Check if (N, d) satisfies the theoretical conditions of Theorem 1 or 2.
    """
    if d % 2 == 0:
        return N >= 3 * d // 2 + 1
    else:
        return N >= 2 * d


def scan_diameter_coverage(n_min=7, n_max=14, tolerance=1e-12):
    """
    Exhaustively check trees for diameter classes 4 <= d <= N - 3.
    
    Returns
    -------
    dict: N -> list of (diameter, has_conflict, is_covered)
    """
    coverage_results = {}
    
    for N in range(n_min, n_max + 1):
        trees_by_diameter = defaultdict(list)
        
        # Enumerate non-isomorphic trees and record (a, W)
        for tree in nx.nonisomorphic_trees(N):
            diam = nx.diameter(tree)
            a_val, w_val = compute_metrics(tree)
            trees_by_diameter[diam].append((a_val, w_val))
            
        row_results = []
        for d in range(4, N - 2):
            trees = trees_by_diameter[d]
            has_conf = has_conflict(trees, tolerance=tolerance)
            covered = is_covered_by_theorems(N, d)
            row_results.append((d, has_conf, covered))
            
        coverage_results[N] = row_results
        
    return coverage_results


def print_coverage_table(results):
    """Print formatted grid of coverage comparison."""
    print(f"{'N':>3} | {'Diameter Coverage [d: Conflict / Theoretical]':<55}")
    print("-" * 65)
    for N, row in results.items():
        formatted_cells = [
            f"d={d}: {'[X]' if c else '[.]'} {'(thm)' if cov else '(-)'}"
            for d, c, cov in row
        ]
        print(f"{N:3d} | {' '.join(formatted_cells)}")


if __name__ == "__main__":
    n_start = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    n_end = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    
    print("=" * 65)
    print("DIAMETER-CLASS CONFLICTS VS THEORETICAL COVERAGE (THEOREMS 1 & 2)")
    print(f"Scanning orders N = {n_start} to {n_end}...")
    print("Legend: [X] = Conflict detected, [.] = No conflict")
    print("        (thm) = Covered by Theorem 1 or 2, (-) = Open/Uncovered")
    print("=" * 65)
    
    results = scan_diameter_coverage(n_start, n_end)
    print_coverage_table(results)
    print("=" * 65)
