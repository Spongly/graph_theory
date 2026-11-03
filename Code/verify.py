r"""
Verification of Theorems 1 and 2 (Conflicts at Every Diameter).

This script tests the explicit conflicting families of trees constructed in:
- Theorem 1 (Even diameter d = 2h, h >= 2):
    A = T(h, h, h, 1^{N - 1 - 3h})
    B = T(h + 1, h - 1, 1^{N - 1 - 2h})
    Claims: a(A) = \lambda_h, a(B) < \lambda_h, W(A) > W(B), diam(A) = diam(B) = 2h.

- Theorem 2 (Odd diameter d = 2h + 1, h >= 2):
    A = T(h + 1, h, h, h, 1^{N - 2 - 4h})
    B = T(h + 2, h - 1, h - 1, 1^{N - 1 - 3h})
    Claims: a(A) >= \lambda_{h+1} > a(B), W(A) > W(B), diam(A) = diam(B) = 2h + 1.

It also verifies the closed-form secular function evaluation \Phi_B(\lambda_h)
against the symbolic formula -\ell * \lambda_h / (1 - \lambda_h).
"""

from math import cos, pi
import sys
from pathlib import Path
import networkx as nx
import numpy as np

# Ensure code directory is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import spider, compute_laplacian_eigenvalues, pinned_path_eigenvalue


def verify_even_diameter_family(h_range=range(2, 10), delta_n=60):
    r"""
    Verify Theorem 1 (even diameter d = 2h).
    Checks that a(A) == \lambda_h, a(B) < \lambda_h, W(A) > W(B), and diam(A) = diam(B) = 2h.
    """
    total_tested = 0
    failures = 0
    
    for h in h_range:
        lam_h = pinned_path_eigenvalue(h)
        for N in range(3 * h + 1, 3 * h + delta_n):
            A = spider([h, h, h], N)
            B = spider([h + 1, h - 1], N)
            
            eig_A = compute_laplacian_eigenvalues(A)
            eig_B = compute_laplacian_eigenvalues(B)
            
            a_A = eig_A[1]
            a_B = eig_B[1]
            w_A = int(nx.wiener_index(A))
            w_B = int(nx.wiener_index(B))
            
            diam_A = nx.diameter(A)
            diam_B = nx.diameter(B)
            
            # Check spectral conditions and Wiener index gap
            is_valid = (
                abs(a_A - lam_h) < 1e-9 and
                a_B < lam_h - 1e-9 and
                w_A > w_B and
                diam_A == 2 * h and
                diam_B == 2 * h
            )
            
            total_tested += 1
            if not is_valid:
                failures += 1
                
    return total_tested, failures


def verify_odd_diameter_family(h_range=range(2, 10), delta_n=60):
    r"""
    Verify Theorem 2 (odd diameter d = 2h + 1).
    Checks that a(A) >= \lambda_{h+1} > a(B), W(A) > W(B), and diam(A) = diam(B) = 2h + 1.
    """
    total_tested = 0
    failures = 0
    
    for h in h_range:
        lam_h1 = pinned_path_eigenvalue(h + 1)
        for N in range(4 * h + 2, 4 * h + delta_n):
            A = spider([h + 1, h, h, h], N)
            B = spider([h + 2, h - 1, h - 1], N)
            
            eig_A = compute_laplacian_eigenvalues(A)
            eig_B = compute_laplacian_eigenvalues(B)
            
            a_A = eig_A[1]
            a_B = eig_B[1]
            w_A = int(nx.wiener_index(A))
            w_B = int(nx.wiener_index(B))
            
            diam_A = nx.diameter(A)
            diam_B = nx.diameter(B)
            
            is_valid = (
                a_A >= lam_h1 - 1e-9 and
                a_B < lam_h1 - 1e-9 and
                w_A > w_B and
                diam_A == 2 * h + 1 and
                diam_B == 2 * h + 1
            )
            
            total_tested += 1
            if not is_valid:
                failures += 1
                
    return total_tested, failures


def evaluate_secular_function(legs, N, theta, lam_val):
    r"""
    Compute secular function \Phi(\lambda) = k - \lambda - \sum g_m(\lambda).
    Here g_m(\lambda) = \cos((m - 0.5)\theta) / \cos((m + 0.5)\theta).
    """
    num_legs = len(legs) + (N - 1 - sum(legs))
    
    def g(m):
        return cos((m - 0.5) * theta) / cos((m + 0.5) * theta)
        
    sum_g = sum(g(m) for m in legs) + (N - 1 - sum(legs)) * g(1)
    return num_legs - lam_val - sum_g


def verify_secular_identities(h_range=range(2, 7)):
    r"""
    Verify that \Phi_B(\lambda_h) matches the closed-form expression:
    -\ell * \lambda_h / (1 - \lambda_h).
    """
    results = []
    for h in h_range:
        theta = pi / (2 * h + 1)
        lam_h = 2.0 - 2.0 * cos(theta)
        N = 3 * h + 6
        ell = N - 1 - 2 * h
        
        phi_numeric = evaluate_secular_function([h + 1, h - 1], N, theta, lam_h)
        phi_closed = -ell * lam_h / (1.0 - lam_h)
        results.append((h, phi_numeric, phi_closed, abs(phi_numeric - phi_closed)))
    return results


if __name__ == "__main__":
    print("=" * 65)
    print("VERIFYING THEOREMS 1 AND 2 ACROSS ORDER N AND DIAMETER d")
    print("=" * 65)
    
    tested_a, failed_a = verify_even_diameter_family()
    print(f"Theorem 1 (Even d): Tested {tested_a} instances, Failures = {failed_a}")
    
    tested_b, failed_b = verify_odd_diameter_family()
    print(f"Theorem 2 (Odd d) : Tested {tested_b} instances, Failures = {failed_b}")
    
    print("\nSecular function Phi_B(lambda_h) comparison with closed form:")
    print(f"{'h':>3} {'Numeric Phi_B':>16} {'Closed Form':>16} {'Diff':>12}")
    print("-" * 50)
    for h, num_val, closed_val, diff in verify_secular_identities():
        print(f"{h:3d} {num_val:16.8f} {closed_val:16.8f} {diff:12.2e}")
    print("=" * 65)
