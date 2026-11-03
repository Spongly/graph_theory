r"""
Generate window_certificates.csv: Exact Rational Certificates for the Window Region.

Proposition (Computer-assisted certificates in the window):
    Let 3 ≤ h ≤ 60, 2h+4 ≤ N ≤ 3h+2, r = N - 2h - 2 (so 2 ≤ r ≤ h), and
        A = T(h+1, h, 2, 1^{r-2})   [legs: [h+1, h, 2] + [1]*(r-2)]
        B = T(h+2, h-1, 1^r)         [legs: [h+2, h-1] + [1]*r]
    Both trees have order N and diameter 2h+1.
    W(A) - W(B) = 2h - 1 - r > 0 (since r ≤ h-1 for the interior range).
    a(B) < a(A) (certified in exact rational arithmetic).

For each (h, r) pair in this range, this script:
1. Builds A and B as starlike (spider) trees.
2. Finds a rational lambda in (a(B), a(A)) using floating-point midpoint.
3. Verifies the certificate via the inertia (Gaussian elimination pivot) lemma.
4. Writes results to data/window_certificates.csv.

Usage:
    python code/generate_window_certificates.py [max_h]
    Default: max_h = 60
"""

import csv
import sys
from fractions import Fraction
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import spider, compute_laplacian_eigenvalues, REPO_ROOT


DATA_DIR = REPO_ROOT / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH = DATA_DIR / "window_certificates.csv"


# ---------------------------------------------------------------------------
# Inertia (pivot) counting for starlike trees
# ---------------------------------------------------------------------------

def spider_pivots(legs: list, N: int, lam: Fraction) -> list:
    """
    Compute all Gaussian-elimination pivots for spider(legs, N) at rational lambda.

    The spider tree has a central vertex (degree = len(legs) + n_pendant) with
    arms of specified lengths, plus n_pendant = N - 1 - sum(legs) pendant leaves.

    Pivot formula (from root outward):
        p_leaf = 1 - lam
        p_v    = deg(v) - lam - sum(1 / p_c  for c children of v)
    """
    one = Fraction(1)
    two = Fraction(2)
    n_pendant = N - 1 - sum(legs)

    pivots = []
    leg_tip_pivots = []   # pivot at the node directly attached to center

    for m in legs:
        p = one - lam   # leaf pivot
        for _ in range(m - 1):
            pivots.append(p)
            p = two - lam - one / p
        leg_tip_pivots.append(p)

    # Pendant leaves (length-1 legs)
    p_pend = one - lam
    for _ in range(n_pendant):
        pivots.append(p_pend)

    # Center pivot
    deg_center = len(legs) + n_pendant
    p_center = Fraction(deg_center) - lam - sum(one / p for p in leg_tip_pivots)
    for p in leg_tip_pivots:
        pivots.append(p)
    pivots.append(p_center)

    return pivots


def count_negative_pivots(legs: list, N: int, lam: Fraction) -> int:
    """Count negative pivots for spider(legs, N) at rational lambda."""
    return sum(1 for p in spider_pivots(legs, N, lam) if p < 0)


def is_valid_certificate(legs_A, legs_B, N, lam):
    """
    Return True if lambda certifies a(A) > a(B):
      - A has exactly 1 negative pivot  => a(A) > lambda
      - B has at least 2 negative pivots => a(B) < lambda
    """
    neg_A = count_negative_pivots(legs_A, N, lam)
    neg_B = count_negative_pivots(legs_B, N, lam)
    return neg_A == 1 and neg_B >= 2


# ---------------------------------------------------------------------------
# Rational separator search
# ---------------------------------------------------------------------------

def find_rational_separator(a_A: float, a_B: float,
                             legs_A: list, legs_B: list, N: int,
                             n_bits: int = 53) -> Fraction:
    """
    Find a rational lambda in (a_B, a_A) that certifies the conflict.

    Start from the floating-point midpoint and verify; if the certificate fails,
    try nearby rationals with denominators that are powers of 2.
    """
    # Try midpoint first
    for bits in range(n_bits, n_bits + 5):
        denom = 2 ** bits
        mid = (a_A + a_B) / 2.0
        num = round(mid * denom)
        lam = Fraction(num, denom)
        if float(lam) > a_B and float(lam) < a_A:
            if is_valid_certificate(legs_A, legs_B, N, lam):
                return lam
    # Fallback: scan numerators near midpoint with larger denominator
    denom = 2 ** (n_bits + 2)
    mid_num = round(((a_A + a_B) / 2.0) * denom)
    for delta in range(0, denom):
        for sign in [0, 1, -1]:
            num = mid_num + sign * delta
            lam = Fraction(num, denom)
            if a_B < float(lam) < a_A:
                if is_valid_certificate(legs_A, legs_B, N, lam):
                    return lam
    raise RuntimeError(f"Could not find rational certificate for N={N}, a_A={a_A}, a_B={a_B}")


# ---------------------------------------------------------------------------
# Main: generate window_certificates.csv
# ---------------------------------------------------------------------------

def generate_window_certificates(max_h: int = 60, max_h_r2: int = 150):
    """
    Generate window_certificates.csv.

    Primary range (Proposition):
        For h = 3..max_h and r = 2..h:
            N = 2h + 2 + r
            A = spider([h+1, h, 2] + [1]*(r-2), N)
            B = spider([h+2, h-1] + [1]*r, N)

    Extended r=2 range (also stated in Proposition):
        For h = max_h+1..max_h_r2 and r = 2 only.
    """
    rows = []
    total = 0
    failed = 0

    # Primary: h=3..max_h, r=2..h
    for h in range(3, max_h + 1):
        for r in range(2, h + 1):
            N = 2 * h + 2 + r
            d = 2 * h + 1

            # Build leg lists
            legs_A = [h + 1, h, 2] + [1] * (r - 2)
            legs_B = [h + 2, h - 1] + [1] * r

            # Verify tree order
            assert 1 + sum(legs_A) + (N - 1 - sum(legs_A)) == N
            assert 1 + sum(legs_B) + (N - 1 - sum(legs_B)) == N

            # Compute a values (floating-point, for initial search)
            A = spider(legs_A, N)
            B = spider(legs_B, N)
            eA = compute_laplacian_eigenvalues(A)
            eB = compute_laplacian_eigenvalues(B)
            a_A = float(eA[1])
            a_B = float(eB[1])

            # Wiener index check
            W_A = int(nx.wiener_index(A))
            W_B = int(nx.wiener_index(B))
            expected_gap = 2 * h - 1 - r
            actual_gap = W_A - W_B

            total += 1

            # Find rational certificate
            try:
                lam = find_rational_separator(a_A, a_B, legs_A, legs_B, N)
                rows.append({
                    "h": h,
                    "r": r,
                    "N": N,
                    "d": d,
                    "lambda_num": lam.numerator,
                    "lambda_den": lam.denominator,
                })
                if h <= 5 or (h % 10 == 0 and r == 2):
                    print(
                        f"  h={h:3d}, r={r:3d}, N={N:4d}, d={d:3d}: "
                        f"a_A={a_A:.8f}, a_B={a_B:.8f}, lam={float(lam):.8f}, "
                        f"W_gap={actual_gap}"
                    )
            except RuntimeError as e:
                print(f"  FAILED h={h}, r={r}: {e}")
                failed += 1

    # Extended r=2 range: h = max_h+1 .. max_h_r2
    for h in range(max_h + 1, max_h_r2 + 1):
        r = 2
        N = 2 * h + 2 + r
        d = 2 * h + 1
        legs_A = [h + 1, h, 2]           # r-2 = 0 extra unit legs
        legs_B = [h + 2, h - 1, 1, 1]   # r = 2 unit legs

        A = spider(legs_A, N)
        B = spider(legs_B, N)
        eA = compute_laplacian_eigenvalues(A)
        eB = compute_laplacian_eigenvalues(B)
        a_A = float(eA[1])
        a_B = float(eB[1])
        W_A = int(nx.wiener_index(A))
        W_B = int(nx.wiener_index(B))

        total += 1
        try:
            lam = find_rational_separator(a_A, a_B, legs_A, legs_B, N)
            rows.append({
                "h": h,
                "r": r,
                "N": N,
                "d": d,
                "lambda_num": lam.numerator,
                "lambda_den": lam.denominator,
            })
            if h % 20 == 1 or h == max_h_r2:
                print(
                    f"  h={h:3d}, r={r:3d}, N={N:4d}, d={d:3d}: "
                    f"a_A={a_A:.8f}, a_B={a_B:.8f}, lam={float(lam):.8f}, "
                    f"W_gap={W_A - W_B}"
                )
        except RuntimeError as e:
            print(f"  FAILED h={h}, r={r}: {e}")
            failed += 1

    print(f"\nGenerated {len(rows)} certificates ({failed} failures out of {total} pairs).")

    with open(OUTPUT_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["h", "r", "N", "d", "lambda_num", "lambda_den"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written: {OUTPUT_PATH} ({len(rows)} rows)")
    return rows


if __name__ == "__main__":
    max_h = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    max_h_r2 = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    print("=" * 70)
    print(f"GENERATING window_certificates.csv (h = 3..{max_h}, r=2..h; plus h={max_h+1}..{max_h_r2} r=2)")
    print("=" * 70)
    generate_window_certificates(max_h, max_h_r2)
    print("=" * 70)
