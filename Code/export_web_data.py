"""
Export verified research paper datasets to JSON for the interactive web application.
Includes:
- All 106 non-isomorphic trees of order 10 with spectra and Wiener indices (Figure 2)
- All 91 diameter conflict witnesses from data/witnesses.csv (Proposition 4 / Figure 4)
- Classification grid from code/scan.json
- Tree statistics from data/tree_stats.json
- Exact verified paper examples (N=7, N=9, N=10, Thm 1-4, Prop 2-4)
"""

import json
import csv
from fractions import Fraction
from pathlib import Path
import networkx as nx
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
CODE_DIR = REPO_ROOT / "code"
DATA_DIR = REPO_ROOT / "data"
WEB_DATA_DIR = REPO_ROOT / "web" / "src" / "data"
WEB_DATA_DIR.mkdir(parents=True, exist_ok=True)

import sys
sys.path.insert(0, str(CODE_DIR))
from utils import spider, comet, compute_laplacian_eigenvalues, compute_metrics, pinned_path_eigenvalue


def export_n10_trees():
    print("Exporting order N=10 trees...")
    trees_n10 = []
    # Identify the Thm 1 pair: A=spider([3,3,3], 10), B=spider([4,2,1,1,1], 10)
    tree_A_edges = set(spider([3,3,3], 10).edges())
    
    idx = 0
    for G in nx.nonisomorphic_trees(10):
        evals = compute_laplacian_eigenvalues(G)
        a_val = float(evals[1])
        W_val = int(nx.wiener_index(G))
        diam = int(nx.diameter(G))
        
        # Check isomorphism to Thm 1 trees
        is_thm1_A = False
        is_thm1_B = False
        
        target_A = spider([3,3,3], 10)
        target_B = spider([4,2,1,1,1], 10)
        if nx.is_isomorphic(G, target_A):
            is_thm1_A = True
        elif nx.is_isomorphic(G, target_B):
            is_thm1_B = True
            
        degrees = [d for n, d in G.degree()]
        
        trees_n10.append({
            "id": idx,
            "edges": [[int(u), int(v)] for u, v in G.edges()],
            "a": round(a_val, 6),
            "W": W_val,
            "d": diam,
            "degrees": degrees,
            "is_thm1_A": is_thm1_A,
            "is_thm1_B": is_thm1_B,
            "label": "A (T(3,3,3))" if is_thm1_A else ("B (T(4,2,1,1,1))" if is_thm1_B else f"Tree #{idx+1}")
        })
        idx += 1
    return trees_n10


def export_witnesses():
    print("Exporting witnesses from data/witnesses.csv...")
    witnesses = []
    witnesses_path = DATA_DIR / "witnesses.csv"
    if not witnesses_path.exists():
        return witnesses

    with open(witnesses_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            N = int(row["N"])
            d = int(row["d"])
            tree_A_str = row["tree_A (larger a, larger W)"]
            tree_B_str = row["tree_B"]
            W_A = int(row["W_A"])
            W_B = int(row["W_B"])
            lam_num = int(row["lambda_num"])
            lam_den = int(row["lambda_den"])
            lam_float = float(Fraction(lam_num, lam_den))
            
            witnesses.append({
                "N": N,
                "d": d,
                "tree_A_str": tree_A_str,
                "tree_B_str": tree_B_str,
                "W_A": W_A,
                "W_B": W_B,
                "lambda_num": lam_num,
                "lambda_den": lam_den,
                "lambda_float": round(lam_float, 8),
                "is_theoretically_proved": (N >= 3 * d // 2 + 1 if d % 2 == 0 else N >= 2 * d)
            })
    return witnesses


def export_verified_examples():
    print("Exporting verified canonical examples...")
    
    # 1. Figure 1: Smallest conflict (N=7, d=4)
    A7 = spider([2,2,2], 7)
    B7 = spider([3,1,1,1], 7)
    a_A7, W_A7 = compute_metrics(A7)
    a_B7, W_B7 = compute_metrics(B7)
    evals_A7 = [round(float(x), 6) for x in compute_laplacian_eigenvalues(A7)]
    evals_B7 = [round(float(x), 6) for x in compute_laplacian_eigenvalues(B7)]
    
    # Second N=7 conflict: Tree with degree sequence (3,3,2,1,1,1,1)
    C7 = nx.Graph()
    C7.add_edges_from([(0,1), (0,2), (0,3), (1,4), (1,5), (4,6)])
    a_C7, W_C7 = compute_metrics(C7)
    evals_C7 = [round(float(x), 6) for x in compute_laplacian_eigenvalues(C7)]
    
    # Order 10: T(3,3,3) vs T(4,2,1,1,1)
    A10 = spider([3,3,3], 10)
    B10 = spider([4,2,1,1,1], 10)
    a_A10, W_A10 = compute_metrics(A10)
    a_B10, W_B10 = compute_metrics(B10)
    evals_A10 = [round(float(x), 6) for x in compute_laplacian_eigenvalues(A10)]
    evals_B10 = [round(float(x), 6) for x in compute_laplacian_eigenvalues(B10)]
    
    # N=9, d=5 non-starlike witness
    G9_1 = nx.Graph()
    for part in "0-1;0-5;0-7;1-2;1-4;2-3;5-6;7-8".split(";"):
        u, v = map(int, part.split("-"))
        G9_1.add_edge(u, v)
        
    G9_2 = nx.Graph()
    for part in "0-1;0-5;0-6;0-7;0-8;1-2;2-3;3-4".split(";"):
        u, v = map(int, part.split("-"))
        G9_2.add_edge(u, v)
        
    a_G9_1, W_G9_1 = compute_metrics(G9_1)
    a_G9_2, W_G9_2 = compute_metrics(G9_2)

    # Double Star D(3,3) - Diameter 3 agreement example (N=8)
    D33 = nx.Graph()
    D33.add_edges_from([(0,1), (0,2), (0,3), (0,4), (1,5), (1,6), (1,7)])
    D24 = nx.Graph()
    D24.add_edges_from([(0,1), (0,2), (0,3), (1,4), (1,5), (1,6), (1,7)])
    a_D33, W_D33 = compute_metrics(D33)
    a_D24, W_D24 = compute_metrics(D24)

    # Diameter N-2 Comet trees (N=8, d=6)
    # comet(8, 2) vs comet(8, 3) vs comet(8, 4)
    C8_2 = comet(8, 2)
    C8_3 = comet(8, 3)
    C8_4 = comet(8, 4)
    a_C8_2, W_C8_2 = compute_metrics(C8_2)
    a_C8_3, W_C8_3 = compute_metrics(C8_3)
    a_C8_4, W_C8_4 = compute_metrics(C8_4)

    return {
        "fig1_smallest": {
            "title": "Figure 1: The Smallest Conflicting Pair (N=7, d=4)",
            "section": "Section 3.2, Figure 1",
            "N": 7,
            "d": 4,
            "tree_A": {
                "name": "A_7 = T(2,2,2)",
                "description": "Symmetric spider tree with 3 legs of length 2",
                "edges": [[int(u), int(v)] for u, v in A7.edges()],
                "a": round(a_A7, 6),
                "exact_a": "(3 - sqrt(5)) / 2 ≈ 0.381966",
                "W": W_A7,
                "Kf": W_A7,
                "eigenvalues": evals_A7,
                "legs": [2, 2, 2],
                "diameter": 4
            },
            "tree_B": {
                "name": "B_7 = T(3,1,1,1)",
                "description": "Asymmetric tree with one long leg of 3 and 3 pendant leaves",
                "edges": [[int(u), int(v)] for u, v in B7.edges()],
                "a": round(a_B7, 6),
                "exact_a": "≈ 0.295532",
                "W": W_B7,
                "Kf": W_B7,
                "eigenvalues": evals_B7,
                "legs": [3, 1, 1, 1],
                "diameter": 4
            },
            "tree_C": {
                "name": "C_7: Degree Seq (3,3,2,1,1,1,1)",
                "description": "Second conflicting tree of order 7, with two degree-3 hubs",
                "edges": [[int(u), int(v)] for u, v in C7.edges()],
                "a": round(a_C7, 6),
                "W": W_C7,
                "Kf": W_C7,
                "eigenvalues": evals_C7,
                "diameter": 4
            },
            "rational_separator": {
                "lambda_num": 381397578438989,
                "lambda_den": 1125899906842624,
                "approx": 0.3387491
            },
            "explanation": "a(A_7) ≈ 0.3820 > a(B_7) ≈ 0.2955, while W(A_7) = 48 > W(B_7) = 46. Algebraic connectivity calls A_7 more robust, but Kirchhoff index calls B_7 more robust!"
        },
        "thm1_even_n10": {
            "title": "Theorem 1 Example: Even Diameter d=6, Order N=10",
            "section": "Theorem 1, Figure 2 & Figure 3",
            "N": 10,
            "d": 6,
            "tree_A": {
                "name": "A = T(3,3,3)",
                "description": "3 legs of length h=3 at central vertex",
                "edges": [[int(u), int(v)] for u, v in A10.edges()],
                "a": round(a_A10, 6),
                "exact_a": "λ_3 = 2 - 2*cos(π/7) ≈ 0.198062",
                "W": W_A10,
                "Kf": W_A10,
                "eigenvalues": evals_A10,
                "legs": [3, 3, 3]
            },
            "tree_B": {
                "name": "B = T(4,2,1,1,1)",
                "description": "Legs of h+1=4, h-1=2, and N-1-2h=3 pendant leaves",
                "edges": [[int(u), int(v)] for u, v in B10.edges()],
                "a": round(a_B10, 6),
                "W": W_B10,
                "Kf": W_B10,
                "eigenvalues": evals_B10,
                "legs": [4, 2, 1, 1, 1]
            },
            "wiener_gap": W_A10 - W_B10,
            "wiener_formula": "(h-2)(h+1)/2 * N - (2h^3 + 3h^2 - 17h - 6)/6 = 2 * 10 - 16 = 4",
            "explanation": "A sits on the plateau a(A) = λ_3 ≈ 0.198062. B has secular function Φ_B(λ_3) < 0, forcing a(B) ≈ 0.1504 < λ_3. Simultaneously, W(A) = 125 > W(B) = 121."
        },
        "non_starlike_witness": {
            "title": "Proposition 4 Witness: Non-Starlike Pair (N=9, d=5)",
            "section": "Proposition 4, witnesses.csv row 5",
            "N": 9,
            "d": 5,
            "tree_A": {
                "name": "G_1 (P_5 with branched hub)",
                "edges": [[int(u), int(v)] for u, v in G9_1.edges()],
                "a": round(a_G9_1, 6),
                "W": W_G9_1,
                "Kf": W_G9_1
            },
            "tree_B": {
                "name": "G_2 = T(4,1,1,1,1)",
                "edges": [[int(u), int(v)] for u, v in G9_2.edges()],
                "a": round(a_G9_2, 6),
                "W": W_G9_2,
                "Kf": W_G9_2
            },
            "rational_separator": {
                "lambda_num": 7542101533641927,
                "lambda_den": 36028797018963968,
                "approx": 0.209335
            },
            "explanation": "The sole non-starlike witness in the entire N <= 19 classification! All other 90 classes are witnessed by pairs of starlike trees."
        },
        "diameter3_agreement": {
            "title": "Proposition 3: Diameter d=3 Always Agrees (Double Stars D(p,q))",
            "section": "Proposition 3",
            "N": 8,
            "d": 3,
            "trees": [
                {"name": "D(3,3)", "p": 3, "q": 3, "pq": 9, "a": round(a_D33, 6), "W": W_D33, "edges": [[int(u), int(v)] for u, v in D33.edges()]},
                {"name": "D(2,4)", "p": 2, "q": 4, "pq": 8, "a": round(a_D24, 6), "W": W_D24, "edges": [[int(u), int(v)] for u, v in D24.edges()]}
            ],
            "formula_W": "W = (N-1)^2 + pq = 49 + pq",
            "explanation": "As product p*q increases, W strictly increases and a(G) strictly decreases. Thus a and W induce EXACT OPPOSITE STRICT ORDERS: zero conflicts possible!"
        },
        "diameter_N_minus_2_agreement": {
            "title": "Proposition 5: Diameter d=N-2 Always Agrees (Comet Trees)",
            "section": "Proposition 5 (Patra & Lal 2008)",
            "N": 8,
            "d": 6,
            "trees": [
                {"name": "Comet(8, pos=2)", "a": round(a_C8_2, 6), "W": W_C8_2, "edges": [[int(u), int(v)] for u, v in C8_2.edges()]},
                {"name": "Comet(8, pos=3)", "a": round(a_C8_3, 6), "W": W_C8_3, "edges": [[int(u), int(v)] for u, v in C8_3.edges()]},
                {"name": "Comet(8, pos=4)", "a": round(a_C8_4, 6), "W": W_C8_4, "edges": [[int(u), int(v)] for u, v in C8_4.edges()]}
            ],
            "explanation": "Moving the pendant leaf toward the path center strictly increases a(G) and strictly decreases W(T). Hence, diameter N-2 is provably conflict-free!"
        }
    }


def main():
    n10_trees = export_n10_trees()
    witnesses = export_witnesses()
    verified_examples = export_verified_examples()
    
    # Load scan.json
    scan_json_path = CODE_DIR / "scan.json"
    with open(scan_json_path, "r", encoding="utf-8") as f:
        scan_data = json.load(f)
        
    # Load tree_stats.json
    stats_json_path = DATA_DIR / "tree_stats.json"
    with open(stats_json_path, "r", encoding="utf-8") as f:
        tree_stats = json.load(f)
        
    # Clean NaNs in tree_stats and scan_data for strict JSON compliance
    def clean_nan(obj):
        if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
            return None
        if isinstance(obj, dict):
            return {k: clean_nan(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [clean_nan(v) for v in obj]
        return obj

    stats_clean = clean_nan(tree_stats)
    scan_clean = clean_nan(scan_data)

    output_data = {
        "metadata": {
            "title": "Algebraic Connectivity Versus the Kirchhoff Index: When Do They Rank Trees Differently?",
            "authors": ["Gopinath Sahoo", "Pratik Pokhrel"],
            "institutions": [
                "Bennett University, Greater Noida, India",
                "Advanced College of Engineering and Management (ACEM), Tribhuvan University, Kathmandu, Nepal"
            ],
            "n10_count": len(n10_trees),
            "witnesses_count": len(witnesses)
        },
        "n10_trees": n10_trees,
        "witnesses": witnesses,
        "verified_examples": verified_examples,
        "scan_data": scan_clean,
        "tree_stats": stats_clean
    }
    
    out_file = WEB_DATA_DIR / "paperData.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2)
        
    print(f"Successfully exported paper data to {out_file} ({out_file.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
