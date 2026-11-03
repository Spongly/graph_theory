r"""
Exhaustive Computations for Small Connected Graphs and Trees (Tables 1 and 2).

This script generates:
- Table 1: Conflicts between algebraic connectivity a(G) and the Kirchhoff index Kf(G)
  among all connected graphs with n \le 7 vertices (from NetworkX graph atlas).
  Includes overall pairs and equal-edge-count pairs.

- Table 2: Conflicting pairs of non-isomorphic trees of order 7 \le N \le 14.
  Uses chunked vectorized matrix broadcasting to count conflicting pairs efficiently.
"""

from pathlib import Path
import sys
import networkx as nx
import numpy as np
from networkx.generators.atlas import graph_atlas_g

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import compute_laplacian_eigenvalues


def compute_graph_spectrum(G: nx.Graph, tolerance=1e-12):
    r"""
    Compute algebraic connectivity a(G), Kirchhoff index Kf(G), and edge count |E|.
    
    By Gutman & Mohar (1996):
        Kf(G) = n * \sum_{i >= 2} 1 / \mu_i
    """
    mu = compute_laplacian_eigenvalues(G)
    n = len(mu)
    a = float(mu[1])
    kf = float(n * np.sum(1.0 / mu[1:]))
    m = G.number_of_edges()
    return a, kf, m


def count_graph_conflicts(graph_specs, tolerance=1e-12):
    """
    Count conflicting pairs among graphs:
    A pair (G_i, G_j) conflicts if (a_i > a_j and Kf_i > Kf_j) or vice versa.
    
    Returns
    -------
    tuple
        (conflicts, total_pairs, equal_edge_conflicts, equal_edge_pairs)
    """
    n_graphs = len(graph_specs)
    if n_graphs < 2:
        return 0, 0, 0, 0
        
    a = np.array([s[0] for s in graph_specs])
    kf = np.array([s[1] for s in graph_specs])
    edges = np.array([s[2] for s in graph_specs])
    
    # Directed differences
    diff_a = a[:, None] - a[None, :]
    diff_kf = kf[:, None] - kf[None, :]
    
    # Conflict matrix: i has strictly larger a and strictly larger Kf
    is_conflict_directed = (diff_a > tolerance) & (diff_kf > tolerance)
    same_edges = (edges[:, None] == edges[None, :])
    
    # Upper triangular mask for unordered pairs
    upper_tri = np.triu(np.ones((n_graphs, n_graphs), dtype=bool), k=1)
    is_conflict_unordered = (is_conflict_directed | is_conflict_directed.T) & upper_tri
    
    total_pairs = n_graphs * (n_graphs - 1) // 2
    total_conflicts = int(is_conflict_unordered.sum())
    
    equal_e_pairs = int((same_edges & upper_tri).sum())
    equal_e_conflicts = int((is_conflict_unordered & same_edges).sum())
    
    return total_conflicts, total_pairs, equal_e_conflicts, equal_e_pairs


def generate_table_1():
    """Generate Table 1 data: Connected graphs of order 4 <= n <= 7."""
    table_rows = []
    for n in range(4, 8):
        connected_graphs = [
            g for g in graph_atlas_g()
            if g.number_of_nodes() == n and nx.is_connected(g)
        ]
        specs = [compute_graph_spectrum(g) for g in connected_graphs]
        conf, pairs, conf_eq, pairs_eq = count_graph_conflicts(specs)
        
        pct = 100.0 * conf / pairs if pairs > 0 else 0.0
        pct_eq = 100.0 * conf_eq / pairs_eq if pairs_eq > 0 else 0.0
        
        table_rows.append({
            "n": n,
            "graphs": len(connected_graphs),
            "pairs": pairs,
            "conflicts": conf,
            "pct": pct,
            "conf_eq": conf_eq,
            "pairs_eq": pairs_eq,
            "pct_eq": pct_eq,
        })
    return table_rows


def generate_table_2(min_N=7, max_N=14, chunk_size=500, tolerance=1e-12):
    """Generate Table 2 data: Conflicting pairs of trees of order min_N <= N <= max_N."""
    table_rows = []
    for N in range(min_N, max_N + 1):
        trees = list(nx.nonisomorphic_trees(N))
        num_trees = len(trees)
        
        # Compute spectral properties
        specs = [compute_graph_spectrum(t) for t in trees]
        a_arr = np.array([s[0] for s in specs])
        w_arr = np.array([round(s[1]) for s in specs], dtype=np.int64)
        
        # Count conflicting pairs (chunked to prevent excessive memory usage)
        total_conflicts = 0
        for i in range(0, num_trees, chunk_size):
            chunk_a = a_arr[i:i + chunk_size, None]
            chunk_w = w_arr[i:i + chunk_size, None]
            
            diff_a = chunk_a - a_arr[None, :]
            diff_w = chunk_w - w_arr[None, :]
            
            is_conf = (diff_a > tolerance) & (diff_w > 0)
            total_conflicts += int(is_conf.sum())
            
        total_pairs = num_trees * (num_trees - 1) // 2
        pct = 100.0 * total_conflicts / total_pairs if total_pairs > 0 else 0.0
        
        table_rows.append({
            "N": N,
            "trees": num_trees,
            "pairs": total_pairs,
            "conflicts": total_conflicts,
            "pct": pct,
        })
    return table_rows


def print_table_1(rows):
    """Print Table 1 in formatted ASCII."""
    print("TABLE 1: Conflicts between a and Kf among connected graphs (n <= 7)")
    print("-" * 75)
    print(f"{'n':>3} {'Graphs':>8} {'Pairs':>10} {'Conflicts':>10} {'%':>8} {'Equal-|E| Conf':>16} {'%':>8}")
    print("-" * 75)
    for r in rows:
        eq_str = f"{r['conf_eq']} / {r['pairs_eq']}"
        print(f"{r['n']:3d} {r['graphs']:8,d} {r['pairs']:10,d} {r['conflicts']:10,d} {r['pct']:8.3f} {eq_str:>16} {r['pct_eq']:8.3f}")
    print("-" * 75)


def print_table_2(rows):
    """Print Table 2 in formatted ASCII."""
    print("\nTABLE 2: Conflicting pairs of non-isomorphic trees (7 <= N <= 14)")
    print("-" * 65)
    print(f"{'N':>3} {'Trees':>8} {'Pairs':>14} {'Conflicts':>14} {'%':>8}")
    print("-" * 65)
    for r in rows:
        print(f"{r['N']:3d} {r['trees']:8,d} {r['pairs']:14,d} {r['conflicts']:14,d} {r['pct']:8.3f}")
    print("-" * 65)


if __name__ == "__main__":
    print("=" * 75)
    print("REPRODUCING TABLES 1 AND 2")
    print("=" * 75)
    
    t1_rows = generate_table_1()
    print_table_1(t1_rows)
    
    # Run Table 2 up to N=12 by default for speed, or sys.argv[1] if specified
    max_tree_n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print(f"\nComputing Tree table up to N = {max_tree_n} (pass argument e.g. '14' for full range):")
    t2_rows = generate_table_2(min_N=7, max_N=max_tree_n)
    print_table_2(t2_rows)
    print("=" * 75)
