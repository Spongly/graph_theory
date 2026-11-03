r"""
Investigation of Diameter N-2 Trees (Comet Trees / Open Problem 1).

In trees of diameter N - 2, vertices form a backbone path P_{N-1} with a single
pendant vertex attached at position i (2 <= i <= floor(N/2)).

Young et al. (2011) proved that the Wiener index W decreases strictly as the
pendant moves toward the center of the path. This script numerically verifies that:
1. Algebraic connectivity a(G) strictly increases toward the center (\Delta a > 0).
2. Wiener index W(G) strictly decreases toward the center (\Delta W < 0).
3. The minimum consecutive gap in a(G) is tracked for N = 7, ..., 120 (approx 8e-9).
"""

from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import comet, compute_metrics


def check_comet_orderings(n_min=7, n_max=120):
    """
    Check whether moving the pendant leaf toward the center strictly increases a
    and strictly decreases W for each order N in range [n_min, n_max].
    
    Returns
    -------
    tuple
        (all_passed: bool, min_a_gap: float, results_summary: list)
    """
    min_consecutive_a_gap = float("inf")
    all_passed = True
    failures = []
    
    for N in range(n_min, n_max + 1):
        backbone_len = N - 1
        # Distinct attachment positions up to symmetry: i = 2 .. floor(N / 2)
        indices = list(range(2, (backbone_len + 1) // 2 + 1))
        
        a_values = []
        w_values = []
        
        for i in indices:
            tree = comet(N, i)
            a_val, w_val = compute_metrics(tree)
            a_values.append(a_val)
            w_values.append(w_val)
            
        diff_a = np.diff(a_values)
        diff_w = np.diff(w_values)
        
        # Verify strict monotonicity toward center
        is_strictly_increasing_a = bool((diff_a > 0).all()) if len(diff_a) > 0 else True
        is_strictly_decreasing_w = bool((diff_w < 0).all()) if len(diff_w) > 0 else True
        
        if len(diff_a) > 0:
            min_consecutive_a_gap = min(min_consecutive_a_gap, float(diff_a.min()))
            
        if not (is_strictly_increasing_a and is_strictly_decreasing_w):
            all_passed = False
            failures.append((N, a_values, w_values))
            
    return all_passed, min_consecutive_a_gap, failures


if __name__ == "__main__":
    print("=" * 65)
    print("VERIFYING COMET TREES (DIAMETER d = N - 2)")
    print("=" * 65)
    print("Checking orders N = 7 to 120...")
    
    success, min_gap, fails = check_comet_orderings(7, 120)
    
    if success:
        print("SUCCESS: For all N in [7, 120], a(T) strictly increases and")
        print("         W(T) strictly decreases toward the center.")
        print(f"Minimum consecutive a-gap observed: {min_gap:.4e}")
    else:
        print(f"FAILURES detected at orders: {[f[0] for f in fails]}")
    print("=" * 65)
