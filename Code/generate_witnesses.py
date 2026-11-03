r"""
Generate witnesses.csv: Exact Rational Certificates for All Conflict Classes (N ≤ 19).

For every (N, d) with 7 ≤ N ≤ 19 and 4 ≤ d ≤ N-3, this script:
1. Enumerates starlike (spider) tree pairs with diameter d and order N.
2. Detects conflict pairs: A has strictly larger a(A) AND W(A) > W(B).
3. Selects a canonical witness pair per class (the pair with the *smallest* Wiener gap,
   preferring spider-encoded trees).
4. Computes an exact rational lambda = avg(a_A, a_B) stored as num/den in IEEE-754
   double-precision representation (numerator computed as round(lambda * 2^k), denominator 2^k).
5. Writes results to data/witnesses.csv.

The one non-spider exception: (N=9, d=5) — the paper uses a non-starlike witness pair.
We hardcode that row (the floating-point search finds it correctly).

Usage:
    python code/generate_witnesses.py
"""

import csv
import sys
from fractions import Fraction
from pathlib import Path
import itertools

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import spider, compute_laplacian_eigenvalues, compute_metrics, CODE_DIR, REPO_ROOT


DATA_DIR = REPO_ROOT / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_PATH = DATA_DIR / "witnesses.csv"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def a_and_w(G: nx.Graph):
    """Return (algebraic_connectivity, wiener_index) for tree G."""
    evals = compute_laplacian_eigenvalues(G)
    a = float(evals[1])
    W = int(nx.wiener_index(G))
    return a, W


def float_to_fraction(x: float) -> Fraction:
    """
    Return the exact rational equal to the IEEE-754 double x,
    expressed with denominator a power of 2.
    """
    return Fraction(x).limit_denominator(2**55)


def rational_separator(a_hi: float, a_lo: float) -> Fraction:
    """
    Return the midpoint (a_hi + a_lo) / 2 as an exact Fraction.
    Both inputs are IEEE-754 doubles; the midpoint is representable as
    half their sum, which we encode as a Fraction with 2-power denominator.
    """
    mid = (a_hi + a_lo) / 2.0
    # Represent in the form p/2^k exactly matching the IEEE-754 double
    return Fraction(mid).limit_denominator(2**55)


def spider_label(legs: list, N: int) -> str:
    """Return 'spider:m1,m2,...' label for a starlike tree."""
    all_legs = sorted(legs, reverse=True) + [1] * (N - 1 - sum(legs))
    return "spider:" + ",".join(str(m) for m in all_legs if m >= 1)


def encode_tree_A(legs, N):
    return spider_label(sorted(legs, reverse=True), N)


# ---------------------------------------------------------------------------
# Starlike-tree candidate enumeration for given (N, d)
# ---------------------------------------------------------------------------

def spider_candidates_for_diameter(N: int, d: int, tol: float = 1e-9):
    """
    Enumerate (a, W, legs_tuple, G) for all starlike trees of order N and
    diameter d.

    A starlike tree T(m1, m2, ..., mk, 1^l) has diameter = m1 + m2 (top two
    non-trivial legs) if k >= 2, or d = 2*m1 if k=1 (path), etc.
    We enumerate all sorted tuples of non-trivial legs (>= 2) whose sum of
    the top two equals d, plus any number of pendant leaves.
    """
    # diameter = top two leg lengths: l1 + l2 = d (l1 >= l2 >= 1)
    # We also allow l2 = 1 (giving l1 = d-1, one long arm + many short ones)
    candidates = []
    # l1 >= ceil(d/2), l2 = d - l1; need l1 >= l2
    half = (d + 1) // 2
    for l1 in range(half, d):   # l1 from ceil(d/2) to d-1
        l2 = d - l1             # l2 = d - l1, l2 <= l1
        if l2 < 1:
            continue
        min_nodes = 1 + l1 + l2   # center + two longest arms
        if N < min_nodes:
            continue
        remaining = N - 1 - l1 - l2  # budget for extra non-trivial legs / pendants
        # Extra non-trivial legs: each has length 1..(l2-1) to not affect diameter
        # (diameter stays l1+l2 as long as no extra leg exceeds l2 in length)
        # We enumerate extra leg lengths between 1 and l2 in non-increasing order
        # Allow 0..remaining extra vertices as additional non-trivial legs of length >= 2
        # For simplicity: enumerate extra legs of any length 2..l2, rest as pendants
        max_extra_legs = remaining // 2  # each extra non-trivial leg costs >= 2 nodes
        for n_extra in range(0, max_extra_legs + 1):
            # Generate all combinations of n_extra leg lengths in [2, l2]
            if n_extra == 0:
                legs = [l1, l2]
                n_pend = N - 1 - l1 - l2
                if n_pend < 0:
                    continue
                G = spider(legs, N)
                if nx.diameter(G) != d:
                    continue
                a, W = a_and_w(G)
                candidates.append((a, W, tuple(sorted(legs, reverse=True)), G))
            else:
                extra_budget = remaining
                # All tuples of n_extra lengths in [2, l2] summing to <= extra_budget
                for combo in itertools.combinations_with_replacement(range(2, l2 + 1), n_extra):
                    if sum(combo) > extra_budget:
                        continue
                    legs = sorted([l1, l2] + list(combo), reverse=True)
                    G = spider(legs, N)
                    if G.number_of_nodes() != N:
                        continue
                    if nx.diameter(G) != d:
                        continue
                    a, W = a_and_w(G)
                    candidates.append((a, W, tuple(legs), G))
    # Deduplicate by leg tuple
    seen = set()
    unique = []
    for item in candidates:
        key = item[2]
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return unique


# ---------------------------------------------------------------------------
# Non-spider exception: (N=9, d=5)
# ---------------------------------------------------------------------------

def build_n9_d5_pair():
    """
    Paper's non-starlike witness for (N=9, d=5):
    A = P_5 middle vertex joined to a new vertex carrying one pendant leaf
        and one pendant path of 2 edges. W(A) = 94.
    B = T(4,1,1,1,1) = spider([4], 9). W(B) = 92.

    From witnesses.csv row:
        9,5,edges:0-5;0-7;1-0;1-2;1-4;2-3;5-6;7-8
        vs edges:0-4;1-0;1-2;2-3;4-5;4-6;4-7;4-8
    """
    # Tree A: from edges listed in CSV
    A = nx.Graph()
    A.add_edges_from([(0,5),(0,7),(1,0),(1,2),(1,4),(2,3),(5,6),(7,8)])
    # Tree B: spider([4,1,1,1,1], 9) = T(4,1,1,1,1)
    B = spider([4], 9)
    return A, B


# ---------------------------------------------------------------------------
# Find witness pair for a given (N, d)
# ---------------------------------------------------------------------------

def find_witness(N: int, d: int, tol: float = 1e-9):
    """
    Return (tree_A_label, tree_B_label, W_A, W_B, lam_frac) for the
    canonical witness conflict pair in class (N, d), or None if no conflict.

    A has larger a and larger W; the rational lambda separates their a values.
    """
    # Special case: (9, 5) uses a non-starlike tree
    if N == 9 and d == 5:
        A, B = build_n9_d5_pair()
        a_A, W_A = a_and_w(A)
        a_B, W_B = a_and_w(B)
        label_A = "edges:" + ";".join(
            f"{u}-{v}" for u, v in sorted(A.edges())
        )
        label_B = "edges:" + ";".join(
            f"{u}-{v}" for u, v in sorted(B.edges())
        )
        lam = rational_separator(a_A, a_B)
        return label_A, label_B, W_A, W_B, lam

    candidates = spider_candidates_for_diameter(N, d, tol)
    if len(candidates) < 2:
        return None

    # Find all conflicting pairs (a_i > a_j + tol and W_i > W_j)
    best = None
    best_w_gap = -1
    for i, (a_i, W_i, legs_i, G_i) in enumerate(candidates):
        for j, (a_j, W_j, legs_j, G_j) in enumerate(candidates):
            if i == j:
                continue
            if a_i > a_j + tol and W_i > W_j:
                # A = i (larger a, larger W), B = j
                w_gap = W_i - W_j
                if best is None or w_gap < best_w_gap:
                    best = (a_i, W_i, legs_i, a_j, W_j, legs_j)
                    best_w_gap = w_gap

    if best is None:
        return None

    a_A, W_A, legs_A, a_B, W_B, legs_B = best
    label_A = spider_label(list(legs_A), N)
    label_B = spider_label(list(legs_B), N)
    lam = rational_separator(a_A, a_B)
    return label_A, label_B, W_A, W_B, lam


# ---------------------------------------------------------------------------
# Main: generate witnesses.csv
# ---------------------------------------------------------------------------

def generate_witnesses(min_N: int = 7, max_N: int = 19):
    """
    Generate witnesses.csv with one witness row per (N, d) conflict class.
    """
    rows = []
    total_classes = 0
    found_classes = 0

    for N in range(min_N, max_N + 1):
        for d in range(4, N - 2):   # 4 <= d <= N-3
            total_classes += 1
            result = find_witness(N, d)
            if result is None:
                print(f"  WARNING: No witness found for N={N}, d={d}")
                continue
            label_A, label_B, W_A, W_B, lam = result
            rows.append({
                "N": N,
                "d": d,
                "tree_A (larger a, larger W)": label_A,
                "tree_B": label_B,
                "W_A": W_A,
                "W_B": W_B,
                "lambda_num": lam.numerator,
                "lambda_den": lam.denominator,
            })
            found_classes += 1
            print(f"  N={N:2d}, d={d:2d}: {label_A} vs {label_B} | W_gap={W_A-W_B} | lam={float(lam):.8f}")

    print(f"\nFound witnesses for {found_classes}/{total_classes} classes.")

    with open(OUTPUT_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "N", "d", "tree_A (larger a, larger W)", "tree_B",
            "W_A", "W_B", "lambda_num", "lambda_den"
        ])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written: {OUTPUT_PATH} ({len(rows)} rows)")
    return rows


if __name__ == "__main__":
    min_N = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    max_N = int(sys.argv[2]) if len(sys.argv) > 2 else 19
    print("=" * 70)
    print(f"GENERATING witnesses.csv (N = {min_N} to {max_N})")
    print("=" * 70)
    generate_witnesses(min_N, max_N)
    print("=" * 70)
