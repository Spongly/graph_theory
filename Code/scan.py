r"""
Exhaustive Diameter-Class Scanner for Trees (Orders N \le 19).

This script performs the diameter-by-diameter classification described in Section 4:
For each order N and diameter d:
1. Generates all non-isomorphic trees of order N.
2. Computes algebraic connectivity a(T) and Wiener index W(T).
3. Evaluates whether the diameter class contains a conflict in O(n \log n) time
   by scanning trees by increasing W and tracking the running minimum of a.
4. Records extremal properties (a_max, W_min) and saves/updates 'scan.json'.

Usage:
    python scan.py <min_N> <max_N>
Example:
    python scan.py 7 10
"""

from collections import defaultdict
import json
from pathlib import Path
import sys
import time
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import CODE_DIR, compute_metrics, has_conflict


def scan_order(N: int, tolerance: float = 1e-12):
    """
    Scan all non-isomorphic trees of order N, grouped by diameter.
    
    Returns
    -------
    dict
        Diameter class dictionary with conflict status and summary statistics.
    """
    trees_by_diameter = defaultdict(list)
    
    # Enumerate all trees of order N
    for tree in nx.nonisomorphic_trees(N):
        diam = nx.diameter(tree)
        a_val, w_val = compute_metrics(tree)
        trees_by_diameter[diam].append((a_val, w_val))
        
    class_summaries = {}
    for diam, items in sorted(trees_by_diameter.items()):
        max_a_item = max(items, key=lambda t: t[0])
        min_w = min(t[1] for t in items)
        
        is_conflicting = has_conflict(items, tolerance=tolerance)
        count_min_w = sum(1 for t in items if t[1] == min_w)
        count_max_a = sum(1 for t in items if abs(t[0] - max_a_item[0]) < tolerance)
        
        class_summaries[diam] = {
            "n": len(items),
            "conflict": is_conflicting,
            "amax_W": max_a_item[1],
            "wmin": min_w,
            "nwmin": count_min_w,
            "namax": count_max_a,
        }
        
    return class_summaries


def run_scanner(min_N: int, max_N: int, json_path: Path = None):
    """
    Run scanner for orders min_N through max_N and save results to JSON.
    """
    if json_path is None:
        json_path = CODE_DIR / "scan.json"
        
    data = {}
    if json_path.exists():
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception as e:
            print(f"Note: Could not load existing {json_path.name} ({e}), starting fresh.")
            
    print(f"Scanning orders N = {min_N} to {max_N}...")
    for N in range(min_N, max_N + 1):
        start_time = time.time()
        order_results = scan_order(N)
        elapsed = time.time() - start_time
        
        # Store using string key for JSON compatibility
        data[str(N)] = order_results
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
            
        conflict_map = {d: int(r["conflict"]) for d, r in sorted(order_results.items())}
        print(f"N = {N:2d} | Elapsed: {elapsed:6.2f}s | Conflicts by diameter: {conflict_map}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        n_start = int(sys.argv[1])
        n_end = int(sys.argv[2])
    elif len(sys.argv) == 2:
        n_start = int(sys.argv[1])
        n_end = n_start
    else:
        n_start, n_end = 7, 8
        print("No order range supplied. Defaulting to N = 7 to 8 (usage: python scan.py <min_N> <max_N>)")
        
    run_scanner(n_start, n_end)
