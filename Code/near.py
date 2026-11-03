r"""
Investigation of Large-Diameter Trees (Diameter d = N - 3 / Open Problem 2).

For trees with diameter close to N, specifically d = N - 3:
- The tree structure consists of a backbone path P_{d+1} (vertices 0..d) plus either:
    1. Two pendant leaves attached at positions i, j along the backbone.
    2. One pendant path of length 2 attached at position i along the backbone.

This script demonstrates that conflicting pairs exist in these large-diameter classes
for orders N = 14, 16, 18, 20, and displays the conflicting pairs with minimal Wiener gaps.
"""

from pathlib import Path
import sys
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import compute_metrics


def build_backbone_tree(diameter: int, branch_spec) -> nx.Graph:
    """
    Construct a tree from a backbone path of length `diameter` (vertices 0..diameter)
    and attach pendant leaves or length-2 paths according to branch_spec.
    
    Parameters
    ----------
    diameter : int
        Number of edges in the backbone path (has diameter + 1 vertices).
    branch_spec : iterable of tuple
        Specifiers: ('leaf', pos) or ('p2', pos), where pos is the backbone attachment node.
    """
    G = nx.path_graph(diameter + 1)
    next_node = diameter + 1
    
    for branch_type, attach_node in branch_spec:
        if branch_type == "leaf":
            G.add_edge(attach_node, next_node)
            next_node += 1
        elif branch_type == "p2":
            G.add_edge(attach_node, next_node)
            G.add_edge(next_node, next_node + 1)
            next_node += 2
            
    return G


def find_diameter_n3_conflicts(orders=(14, 16, 18, 20), tolerance=1e-12):
    """
    Find and return conflicting pairs of trees with diameter d = N - 3.
    """
    results = {}
    
    for N in orders:
        d = N - 3
        # Candidate configurations: two pendant leaves or one pendant path of length 2
        specs = (
            [(("leaf", i), ("leaf", j)) for i in range(1, d) for j in range(i, d)] +
            [(("p2", i),) for i in range(2, d - 1)]
        )
        
        candidates = []
        for spec in specs:
            tree = build_backbone_tree(d, spec)
            if nx.diameter(tree) != d:
                continue
            a_val, w_val = compute_metrics(tree)
            candidates.append((spec, a_val, w_val))
            
        # Detect conflicting pairs: tree_x has strictly larger a AND strictly larger W
        conflicting_pairs = [
            (tree_x, tree_y)
            for tree_x in candidates
            for tree_y in candidates
            if tree_x[1] > tree_y[1] + tolerance and tree_x[2] > tree_y[2]
        ]
        
        # Sort by Wiener gap ascending (smallest gap first)
        conflicting_pairs.sort(key=lambda pair: (pair[0][2] - pair[1][2]))
        results[N] = {
            "diameter": d,
            "num_candidates": len(candidates),
            "num_conflicts": len(conflicting_pairs),
            "top_pairs": conflicting_pairs[:6],
        }
        
    return results


if __name__ == "__main__":
    print("=" * 70)
    print("LARGE DIAMETER ANALYSIS: TREES WITH d = N - 3 (OPEN PROBLEM 2)")
    print("=" * 70)
    
    conflict_data = find_diameter_n3_conflicts()
    
    for N, data in conflict_data.items():
        d = data["diameter"]
        n_trees = data["num_candidates"]
        n_conflicts = data["num_conflicts"]
        
        print(f"\nOrder N = {N:2d} (Diameter d = {d:2d}) | Candidates: {n_trees:3d} | Conflicting pairs: {n_conflicts}")
        print("  Top conflicting pairs with smallest Wiener index gap Delta W:")
        
        for tree_x, tree_y in data["top_pairs"]:
            spec_x, a_x, w_x = tree_x
            spec_y, a_y, w_y = tree_y
            gap_w = w_x - w_y
            gap_a = a_x - a_y
            print(f"   - {str(spec_x):<25} (W={w_x:3d}, a={a_x:.5f}) vs "
                  f"{str(spec_y):<25} (W={w_y:3d}, a={a_y:.5f}) | "
                  f"Delta W = {gap_w:2d}, Delta a = +{gap_a:.5f}")
                  
    print("=" * 70)
