r"""
Verify witnesses.csv and window_certificates.csv using the Inertia Lemma.

Lemma (Exact eigenvalue counting on trees):
    Let T be a tree, rooted at some vertex, and let lambda be a rational number.
    Process the vertices from leaves towards the root and compute pivots:
        p_v = deg(v) - lambda - sum(1/p_c  for each child c of v)
    If no pivot is zero, then lambda is NOT a Laplacian eigenvalue of T, and the
    number of Laplacian eigenvalues of T strictly less than lambda equals the
    number of vertices with p_v < 0.

Certificate interpretation:
    A rational lambda certifies a(A) > a(B) if:
    - A has exactly 1 negative pivot at lambda  =>  exactly one eigenvalue < lambda  =>  a(A) > lambda
    - B has at least 2 negative pivots at lambda =>  at least two eigenvalues < lambda  =>  a(B) < lambda

This script:
1. Reads data/witnesses.csv and verifies each row's rational lambda.
2. Reads data/window_certificates.csv and verifies each row's rational lambda.
3. Prints a full summary with pass/fail counts.

Usage:
    python code/verify_certificates.py

Exit code 0 if all certificates pass; 1 otherwise.
"""

import csv
import sys
from fractions import Fraction
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import spider, REPO_ROOT


DATA_DIR = REPO_ROOT / "data"
WITNESSES_PATH = DATA_DIR / "witnesses.csv"
WINDOW_PATH = DATA_DIR / "window_certificates.csv"


# ---------------------------------------------------------------------------
# Inertia pivot computation (exact rational arithmetic)
# ---------------------------------------------------------------------------

def spider_pivots(legs: list, N: int, lam: Fraction) -> list:
    """
    Compute all Gaussian-elimination pivots for spider(legs, N) at rational lambda.

    Pivot formula (bottom-up from leaves to center):
        p_leaf = 1 - lam
        p_v    = deg(v) - lam - sum(1/p_c  for children c)
    """
    one = Fraction(1)
    two = Fraction(2)
    n_pendant = N - 1 - sum(legs)

    pivots = []
    leg_tip_pivots = []

    for m in legs:
        p = one - lam          # leaf
        for _ in range(m - 1):
            pivots.append(p)
            p = two - lam - one / p
        leg_tip_pivots.append(p)   # node directly attached to center

    p_pend = one - lam
    for _ in range(n_pendant):
        pivots.append(p_pend)

    deg_center = len(legs) + n_pendant
    p_center = Fraction(deg_center) - lam - sum(one / p for p in leg_tip_pivots)
    for p in leg_tip_pivots:
        pivots.append(p)
    pivots.append(p_center)

    return pivots


def general_tree_pivots(G: nx.Graph, lam: Fraction) -> list:
    """
    Compute pivots for a general (non-spider) tree G at rational lambda,
    rooted at an arbitrary vertex via BFS.
    """
    root = next(iter(G.nodes()))
    parent = {root: None}
    order = []   # BFS order
    queue = [root]
    while queue:
        v = queue.pop(0)
        order.append(v)
        for u in G.neighbors(v):
            if u not in parent:
                parent[u] = v
                queue.append(u)

    one = Fraction(1)
    p = {}
    child_sum = {v: Fraction(0) for v in G.nodes()}

    for v in reversed(order):
        deg = Fraction(G.degree(v))
        p[v] = deg - lam - child_sum[v]
        par = parent[v]
        if par is not None:
            child_sum[par] += one / p[v]

    return list(p.values())


def count_negative_pivots_spider(legs: list, N: int, lam: Fraction) -> int:
    return sum(1 for pv in spider_pivots(legs, N, lam) if pv < 0)


def count_negative_pivots_general(G: nx.Graph, lam: Fraction) -> int:
    return sum(1 for pv in general_tree_pivots(G, lam) if pv < 0)


def check_no_zero_pivot(pivots: list) -> bool:
    return all(pv != 0 for pv in pivots)


# ---------------------------------------------------------------------------
# Parse tree specifications from CSV
# ---------------------------------------------------------------------------

def parse_spider_legs(spec: str) -> list:
    """Parse 'spider:a,b,c,...' into leg list [a,b,c,...]."""
    assert spec.startswith("spider:"), f"Expected spider spec, got: {spec}"
    parts = spec[len("spider:"):].split(",")
    legs = [int(x) for x in parts if int(x) >= 1]
    return legs


def parse_edge_tree(spec: str) -> nx.Graph:
    """Parse 'edges:u1-v1;u2-v2;...' into a networkx graph."""
    assert spec.startswith("edges:"), f"Expected edges spec, got: {spec}"
    G = nx.Graph()
    for edge_str in spec[len("edges:"):].split(";"):
        u, v = edge_str.split("-")
        G.add_edge(int(u), int(v))
    return G


# ---------------------------------------------------------------------------
# Verify witnesses.csv
# ---------------------------------------------------------------------------

def verify_witnesses(path: Path = WITNESSES_PATH, verbose: bool = True) -> dict:
    """
    Verify every certificate in witnesses.csv.

    For each row, check:
    1. W_A > W_B  (Wiener index ordering)
    2. lambda strictly between a(B) and a(A) according to pivot counts:
       - neg_pivots(A, lam) == 1
       - neg_pivots(B, lam) >= 2
    """
    if not path.exists():
        print(f"ERROR: {path} not found.")
        return {"passed": 0, "failed": 1, "total": 1}

    passed = 0
    failed = 0
    failures = []

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    for row in rows:
        N = int(row["N"])
        d = int(row["d"])
        spec_A = row["tree_A (larger a, larger W)"].strip()
        spec_B = row["tree_B"].strip()
        W_A = int(row["W_A"])
        W_B = int(row["W_B"])
        lam = Fraction(int(row["lambda_num"]), int(row["lambda_den"]))

        # Build trees and compute pivot counts
        is_spider_A = spec_A.startswith("spider:")
        is_spider_B = spec_B.startswith("spider:")

        ok = True
        reason = ""

        # Check Wiener ordering
        if W_A <= W_B:
            ok = False
            reason = f"W_A={W_A} not > W_B={W_B}"

        if ok:
            # Compute pivot counts
            if is_spider_A:
                legs_A = parse_spider_legs(spec_A)
                pA = spider_pivots(legs_A, N, lam)
                neg_A = sum(1 for pv in pA if pv < 0)
                zero_A = not check_no_zero_pivot(pA)
            else:
                G_A = parse_edge_tree(spec_A)
                pA = general_tree_pivots(G_A, lam)
                neg_A = sum(1 for pv in pA if pv < 0)
                zero_A = not check_no_zero_pivot(pA)

            if is_spider_B:
                legs_B = parse_spider_legs(spec_B)
                pB = spider_pivots(legs_B, N, lam)
                neg_B = sum(1 for pv in pB if pv < 0)
                zero_B = not check_no_zero_pivot(pB)
            else:
                G_B = parse_edge_tree(spec_B)
                pB = general_tree_pivots(G_B, lam)
                neg_B = sum(1 for pv in pB if pv < 0)
                zero_B = not check_no_zero_pivot(pB)

            if zero_A or zero_B:
                ok = False
                reason = f"Zero pivot encountered (lambda is an eigenvalue)"
            elif neg_A != 1:
                ok = False
                reason = f"A has {neg_A} negative pivots (expected 1)"
            elif neg_B < 2:
                ok = False
                reason = f"B has {neg_B} negative pivots (expected >= 2)"

        if ok:
            passed += 1
            if verbose:
                print(f"  PASS N={N:2d}, d={d:2d}: {spec_A} vs {spec_B}")
        else:
            failed += 1
            failures.append((N, d, reason))
            print(f"  FAIL N={N:2d}, d={d:2d}: {reason}")

    return {"passed": passed, "failed": failed, "total": passed + failed, "failures": failures}


# ---------------------------------------------------------------------------
# Verify window_certificates.csv
# ---------------------------------------------------------------------------

def verify_window_certificates(path: Path = WINDOW_PATH, verbose: bool = False) -> dict:
    """
    Verify every certificate in window_certificates.csv.

    For each row (h, r, N, d, lambda_num, lambda_den):
        A = spider([h+1, h, 2, 1*(r-2)], N)
        B = spider([h+2, h-1, 1*r], N)
    Check:
        neg_pivots(A, lam) == 1  and  neg_pivots(B, lam) >= 2
    """
    if not path.exists():
        print(f"ERROR: {path} not found.")
        return {"passed": 0, "failed": 1, "total": 1}

    passed = 0
    failed = 0
    failures = []

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    for row in rows:
        h = int(row["h"])
        r = int(row["r"])
        N = int(row["N"])
        d = int(row["d"])
        lam = Fraction(int(row["lambda_num"]), int(row["lambda_den"]))

        legs_A = [h + 1, h, 2] + [1] * (r - 2)
        legs_B = [h + 2, h - 1] + [1] * r

        pA = spider_pivots(legs_A, N, lam)
        pB = spider_pivots(legs_B, N, lam)

        neg_A = sum(1 for pv in pA if pv < 0)
        neg_B = sum(1 for pv in pB if pv < 0)
        zero_piv = not check_no_zero_pivot(pA) or not check_no_zero_pivot(pB)

        ok = (not zero_piv) and (neg_A == 1) and (neg_B >= 2)
        d_expected = 2 * h + 1
        if d != d_expected:
            ok = False
            reason = f"d={d} != 2h+1={d_expected}"
        elif zero_piv:
            reason = "Zero pivot encountered"
        elif neg_A != 1:
            reason = f"A has {neg_A} neg pivots (expected 1)"
        elif neg_B < 2:
            reason = f"B has {neg_B} neg pivots (expected >= 2)"
        else:
            reason = "OK"

        if ok:
            passed += 1
            if verbose:
                print(f"  PASS h={h:3d}, r={r:3d}, N={N:4d}, d={d:3d}")
        else:
            failed += 1
            failures.append((h, r, N, d, reason))
            print(f"  FAIL h={h:3d}, r={r:3d}, N={N:4d}, d={d:3d}: {reason}")

    return {"passed": passed, "failed": failed, "total": passed + failed, "failures": failures}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 70)
    print("VERIFYING EXACT RATIONAL CERTIFICATES (INERTIA LEMMA)")
    print("=" * 70)

    print(f"\n[1/2] Verifying witnesses.csv ({WITNESSES_PATH.name})")
    print("-" * 60)
    w_result = verify_witnesses(verbose=True)
    print(f"\nResult: {w_result['passed']}/{w_result['total']} passed"
          + (f", {w_result['failed']} FAILURES" if w_result['failed'] else " — ALL OK"))

    print(f"\n[2/2] Verifying window_certificates.csv ({WINDOW_PATH.name})")
    print("-" * 60)
    c_result = verify_window_certificates(verbose=False)
    print(f"Checked {c_result['total']} certificates: "
          f"{c_result['passed']} passed, {c_result['failed']} failed.")
    if c_result['failed'] == 0:
        print("  ALL window certificates PASSED.")
    else:
        print("  FAILURES:")
        for entry in c_result['failures'][:10]:
            print(f"    h={entry[0]}, r={entry[1]}, N={entry[2]}, d={entry[3]}: {entry[4]}")

    print("\n" + "=" * 70)
    all_ok = w_result['failed'] == 0 and c_result['failed'] == 0
    if all_ok:
        print("ALL CERTIFICATES VERIFIED SUCCESSFULLY.")
    else:
        total_fail = w_result['failed'] + c_result['failed']
        print(f"VERIFICATION FAILED: {total_fail} certificate(s) did not pass.")
    print("=" * 70)

    sys.exit(0 if all_ok else 1)
