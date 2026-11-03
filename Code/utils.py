"""
Common utility functions and helpers for tree generation and spectral metrics.

This module provides reusable routines for:
- Constructing starlike (spider) and comet trees
- Computing Laplacian spectra, algebraic connectivity a(G), and Wiener index W(G)
- Standardized directory paths for figures and data files
"""

from pathlib import Path
from typing import List, Tuple
import networkx as nx
import numpy as np

# Standard directories (independent of current working directory)
CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent
FIGURES_DIR = REPO_ROOT / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def spider(legs: List[int], N: int) -> nx.Graph:
    """
    Construct a starlike (spider) tree T(m_1, ..., m_k, 1^l).
    
    Parameters
    ----------
    legs : list of int
        Lengths of non-trivial branches [m_1, ..., m_k] radiating from central vertex 0.
    N : int
        Total number of vertices in the tree. The remaining vertices
        l = N - 1 - sum(legs) are attached as pendant leaves (branches of length 1).
    
    Returns
    -------
    nx.Graph
        The constructed tree on N vertices.
    """
    G = nx.Graph()
    G.add_node(0)
    current_node = 1
    
    for length in legs:
        prev = 0
        for _ in range(length):
            G.add_edge(prev, current_node)
            prev = current_node
            current_node += 1
            
    # Attach remaining vertices as pendant leaves to the center (node 0)
    while current_node < N:
        G.add_edge(0, current_node)
        current_node += 1
        
    return G


def comet(N: int, i: int) -> nx.Graph:
    """
    Construct a comet tree: a path P_{N-1} with a single pendant leaf attached at vertex i.
    
    Parameters
    ----------
    N : int
        Total number of vertices.
    i : int
        1-based position along the backbone path (2 <= i <= (N-1)//2).
        The pendant vertex (index N-1) is attached to node i-1.
        
    Returns
    -------
    nx.Graph
        Comet graph on N vertices.
    """
    G = nx.path_graph(N - 1)
    G.add_edge(i - 1, N - 1)
    return G


def compute_laplacian_eigenvalues(G: nx.Graph) -> np.ndarray:
    r"""
    Compute sorted Laplacian eigenvalues of a graph: 0 = \mu_1 <= \mu_2 <= ... <= \mu_n.
    """
    L = nx.laplacian_matrix(G).toarray().astype(float)
    return np.linalg.eigvalsh(L)


def compute_algebraic_connectivity(G: nx.Graph) -> float:
    r"""
    Compute algebraic connectivity a(G) = \mu_2(G) (the second-smallest Laplacian eigenvalue).
    """
    evals = compute_laplacian_eigenvalues(G)
    return float(evals[1])


def compute_wiener_index(G: nx.Graph) -> int:
    r"""
    Compute the Wiener index W(G) = \sum_{u < v} d(u, v).
    """
    return int(nx.wiener_index(G))


def compute_metrics(G: nx.Graph) -> Tuple[float, int]:
    r"""
    Compute (algebraic_connectivity, wiener_index) for a connected tree.
    
    Uses Gutman & Mohar identity: W(T) = Kf(T) = n * \sum_{i >= 2} 1 / \mu_i.
    """
    evals = compute_laplacian_eigenvalues(G)
    a = float(evals[1])
    n = G.number_of_nodes()
    W = int(round(n * np.sum(1.0 / evals[1:])))
    return a, W


def pinned_path_eigenvalue(m: int) -> float:
    r"""
    Compute smallest eigenvalue \lambda_m = 2 - 2*cos(\pi / (2m + 1)) of pinned path \widehat{P}_m.
    """
    return float(2.0 - 2.0 * np.cos(np.pi / (2.0 * m + 1.0)))


def has_conflict(items: List[Tuple[float, int]], tolerance: float = 1e-12) -> bool:
    r"""
    Determine whether a collection of (a, W) pairs contains a conflicting pair.
    
    A conflict occurs if there exist x, y such that W_x > W_y and a_x > a_y + tol.
    
    Algorithm (O(n log n)):
    1. Sort items by W ascending.
    2. Group items with identical W.
    3. Maintain `running_min_a` as the minimum `a` among all items with strictly smaller W.
    4. If any item in the current group has a > running_min_a + tol, a conflict exists.
    """
    if len(items) < 2:
        return False
        
    sorted_items = sorted(items, key=lambda t: t[1])
    running_min_a = float("inf")
    idx = 0
    total = len(sorted_items)
    
    while idx < total:
        grp_end = idx
        while grp_end < total and sorted_items[grp_end][1] == sorted_items[idx][1]:
            grp_end += 1
            
        group_a = [t[0] for t in sorted_items[idx:grp_end]]
        if max(group_a) > running_min_a + tolerance:
            return True
            
        running_min_a = min(running_min_a, min(group_a))
        idx = grp_end
        
    return False
