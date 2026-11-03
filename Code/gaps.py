r"""
Symbolic derivation and numerical verification of Wiener index gaps.

This script symbolically computes the Wiener index gaps from Theorems 1 and 2:
- Theorem 1:
    W(A) - W(B) = \frac{(h - 2)(h + 1)}{2} N - \frac{2h^3 + 3h^2 - 17h - 6}{6}
- Theorem 2:
    W(A) - W(B) = \frac{(h - 2)(h + 3)}{2} N - \frac{2h^3 + 9h^2 - 29h - 30}{6}

It then asserts that the closed-form symbolic polynomials exactly match the
Wiener indices computed via NetworkX on the actual graph models.
"""

from pathlib import Path
import sys
import sympy as sp
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import spider


def symbolic_omega(m, N):
    r"""
    Wiener contribution of a leg of length m in a tree of order N (Lemma 5):
    \omega_N(m) = N * m*(m+1)/2 - m*(m+1)*(2m+1)/6.
    """
    return N * (m * (m + 1)) / 2 - (m * (m + 1) * (2 * m + 1)) / 6


def symbolic_wiener_starlike(legs, N):
    r"""
    Total Wiener index of starlike tree T(m_1, ..., m_k, 1^\ell).
    Leaves (length 1) have \omega_N(1) = N - 1.
    """
    sum_legs = sum(legs)
    num_leaves = N - 1 - sum_legs
    return sum(symbolic_omega(m, N) for m in legs) + num_leaves * (N - 1)


def derive_symbolic_gaps():
    r"""
    Derive simplified polynomials for \Delta W_A and \Delta W_B.
    """
    h, N = sp.symbols('h N')
    
    # Theorem 1: A = T(h, h, h), B = T(h+1, h-1)
    gap_A = sp.expand(
        symbolic_wiener_starlike([h, h, h], N) -
        symbolic_wiener_starlike([h + 1, h - 1], N)
    )
    
    # Theorem 2: A = T(h+1, h, h, h), B = T(h+2, h-1, h-1)
    gap_B = sp.expand(
        symbolic_wiener_starlike([h + 1, h, h, h], N) -
        symbolic_wiener_starlike([h + 2, h - 1, h - 1], N)
    )
    
    return h, N, gap_A, gap_B


def verify_with_networkx(h_sym, N_sym, gap_A_sym, gap_B_sym, h_vals=range(2, 6)):
    """
    Verify that symbolic polynomial predictions match concrete graph calculations.
    """
    for hh in h_vals:
        for NN in (4 * hh + 2, 4 * hh + 9):
            # Graphs for Theorem 1
            tree_A_even = spider([hh, hh, hh], NN)
            tree_B_even = spider([hh + 1, hh - 1], NN)
            diff_even = int(nx.wiener_index(tree_A_even) - nx.wiener_index(tree_B_even))
            predicted_even = int(gap_A_sym.subs({h_sym: hh, N_sym: NN}))
            
            # Graphs for Theorem 2
            tree_A_odd = spider([hh + 1, hh, hh, hh], NN)
            tree_B_odd = spider([hh + 2, hh - 1, hh - 1], NN)
            diff_odd = int(nx.wiener_index(tree_A_odd) - nx.wiener_index(tree_B_odd))
            predicted_odd = int(gap_B_sym.subs({h_sym: hh, N_sym: NN}))
            
            assert diff_even == predicted_even, f"Even mismatch at h={hh}, N={NN}"
            assert diff_odd == predicted_odd, f"Odd mismatch at h={hh}, N={NN}"


if __name__ == "__main__":
    print("=" * 65)
    print("SYMBOLIC DERIVATION OF WIENER INDEX GAPS")
    print("=" * 65)
    
    h, N, gap_A, gap_B = derive_symbolic_gaps()
    
    print("Theorem 1 Wiener gap  Delta W_A:")
    print("  6 * Delta W_A =", sp.factor(6 * gap_A))
    print("  At threshold N = 3h + 1:", sp.factor(gap_A.subs(N, 3 * h + 1)))
    
    print("\nTheorem 2 Wiener gap  Delta W_B:")
    print("  6 * Delta W_B =", sp.factor(6 * gap_B))
    print("  At threshold N = 4h + 2:", sp.factor(gap_B.subs(N, 4 * h + 2)))
    
    print("\nVerifying symbolic expressions against NetworkX graph computations...")
    verify_with_networkx(h, N, gap_A, gap_B)
    print("PASS: All symbolic gap formulas match NetworkX values exactly.")
    print("=" * 65)
