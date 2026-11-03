"""Generates fig1_smallest_pair.pdf, fig2_scatter_n10.pdf, fig3_AB_curves.pdf (all data computed, nothing hard-coded)
and prints the numbers quoted in the paper (Remark 1, Section 3.2)."""
import numpy as np, networkx as nx, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from exact import spider
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
def a_of(T): return float(np.linalg.eigvalsh(nx.laplacian_matrix(T).toarray().astype(float))[1])
def W_of(T): return int(nx.wiener_index(T))
lam = lambda m: 2 - 2*np.cos(np.pi/(2*m + 1))

# ---------------- Figure 1: the smallest conflicting pair (N = 7) ----------------
def spider_layout(legs, base_angle=90):
    pos = {0: (0.0, 0.0)}; n = 1; k = len(legs)
    for i, m in enumerate(legs):
        ang = np.deg2rad(base_angle + 360.0*i/k)
        for j in range(1, m + 1):
            pos[n] = (j*np.cos(ang), j*np.sin(ang)); n += 1
    return pos
A7, B7 = spider([2, 2, 2]), spider([3, 1, 1, 1])
fig, axs = plt.subplots(1, 2, figsize=(6.0, 2.5))
for ax, T, legs, name in ((axs[0], A7, [2, 2, 2], "A_7"), (axs[1], B7, [3, 1, 1, 1], "B_7")):
    pos = spider_layout(legs, 90 if len(legs) == 3 else 90)
    nx.draw_networkx_edges(T, pos, ax=ax, width=1.6, edge_color="#444444")
    nx.draw_networkx_nodes(T, pos, ax=ax, node_size=110, node_color=["#f28e2b" if v == 0 else "#3b6ea5" for v in T.nodes], edgecolors="k", linewidths=.6)
    ax.set_title(r"$%s=T(%s)$" % (name, ",".join(map(str, legs))), fontsize=10)
    ax.text(0.5, -0.04, r"$a=%.4f$,  $W=%d$" % (a_of(T), W_of(T)), transform=ax.transAxes, ha="center", va="top", fontsize=9)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-3.0, 3.0); ax.set_ylim(-1.5, 3.3)
fig.tight_layout(); fig.savefig("fig1_smallest_pair.pdf"); plt.close(fig)
print("Fig 1: a(A7)=%.4f W=%d | a(B7)=%.4f W=%d | (3-sqrt5)/2=%.4f" % (a_of(A7), W_of(A7), a_of(B7), W_of(B7), (3-5**.5)/2))
# the other conflicting partner of A7 mentioned in Section 3.2
for T in nx.nonisomorphic_trees(7):
    if sorted(dict(T.degree()).values(), reverse=True) == [3, 3, 2, 1, 1, 1, 1]: print("   degree sequence (3,3,2,1,1,1,1): a=%.4f W=%d" % (a_of(T), W_of(T)))

# ---------------- Figure 2: all trees of order 10 ----------------
trees = list(nx.nonisomorphic_trees(10)); a = np.array([a_of(T) for T in trees]); W = np.array([W_of(T) for T in trees]); D = np.array([nx.diameter(T) for T in trees])
n = len(trees); inpair = np.zeros(n, bool)
for i in range(n):
    c = (a[i] > a + 1e-12) & (W[i] > W); inpair[i] |= c.any(); inpair[c] = True
print("Fig 2: trees of order 10:", n, "| trees in at least one conflicting pair:", int(inpair.sum()))
TA, TB = spider([3, 3, 3]), spider([4, 2, 1, 1, 1]); aA, wA, aB, wB = a_of(TA), W_of(TA), a_of(TB), W_of(TB)
print("   A=T(3,3,3): a=%.4f W=%d | B=T(4,2,1,1,1): a=%.4f W=%d | conflict: %s" % (aA, wA, aB, wB, aA > aB and wA > wB))
fig, ax = plt.subplots(figsize=(5.6, 3.9)); cm = plt.get_cmap("viridis")
for d in sorted(set(D)):
    m = D == d; ax.scatter(W[m], a[m], s=16, color=cm((d - 2)/7.0), edgecolors="k", linewidths=.3, label=r"$d=%d$" % d, zorder=3)
for (x, y, lab, dx, dy) in ((wA, aA, r"$A=T(3,3,3)$", 4, 0.035), (wB, aB, r"$B=T(4,2,1,1,1)$", -1.5, -0.075)):
    ax.scatter([x], [y], s=150, facecolors="none", edgecolors="crimson", linewidths=1.5, zorder=4)
    ax.annotate(lab, (x, y), xytext=(x + dx, y + dy), fontsize=8.5, color="crimson", arrowprops=dict(arrowstyle="-", color="crimson", lw=.8))
ax.set_xlabel(r"Wiener index $W(T)$"); ax.set_ylabel(r"algebraic connectivity $a(T)$"); ax.legend(ncol=2, fontsize=7.5, frameon=False, title="diameter", title_fontsize=8)
fig.tight_layout(); fig.savefig("fig2_scatter_n10.pdf"); plt.close(fig)

# ---------------- Figure 3: a(A), a(B) of Theorems 1 and 3 ----------------
cols = {2: "#3b6ea5", 3: "#f28e2b", 4: "#59a14f"}
fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.1))
for h in (2, 3, 4):
    Ns = list(range(3*h + 1, 3*h + 41))
    aa = [a_of(spider([h, h, h] + [1]*(N - 1 - 3*h))) for N in Ns]; ab = [a_of(spider([h + 1, h - 1] + [1]*(N - 1 - 2*h))) for N in Ns]
    assert all(abs(x - lam(h)) < 1e-9 for x in aa) and all(y < lam(h) for y in ab)
    axs[0].plot(Ns, aa, "o", ms=2.6, color=cols[h]); axs[0].plot(Ns, ab, "-", color=cols[h], lw=1.2, label=r"$h=%d$" % h)
    axs[0].axhline(lam(h), color=cols[h], ls="--", lw=.8)
    Ns = list(range(4*h + 2, 4*h + 42))
    aa = [a_of(spider([h + 1, h, h, h] + [1]*(N - 2 - 4*h))) for N in Ns]; ab = [a_of(spider([h + 2, h - 1, h - 1] + [1]*(N - 1 - 3*h))) for N in Ns]
    assert all(x >= lam(h + 1) for x in aa) and all(y < lam(h + 1) for y in ab)
    axs[1].plot(Ns, aa, "o", ms=2.6, color=cols[h]); axs[1].plot(Ns, ab, "-", color=cols[h], lw=1.2, label=r"$h=%d$" % h)
    axs[1].axhline(lam(h + 1), color=cols[h], ls="--", lw=.8)
axs[0].set_title(r"Theorem 1 ($h=2,3,4$)", fontsize=9)
axs[1].set_title(r"Theorem 3 ($h=2,3,4$)", fontsize=9)
from matplotlib.lines import Line2D
for ax in axs: ax.set_xlabel(r"order $N$"); ax.set_ylabel(r"algebraic connectivity")
handles = [Line2D([], [], color=cols[h], lw=2.6, label=r"$h=%d$" % h) for h in (2, 3, 4)] + [
    Line2D([], [], color="k", marker="o", ms=3, ls="", label=r"$a(A)$ (dots)"), Line2D([], [], color="k", lw=1.2, label=r"$a(B)$ (solid)"),
    Line2D([], [], color="k", lw=.8, ls="--", label=r"plateau $\lambda_h$ (left), $\lambda_{h+1}$ (right)")]
fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=7.5, frameon=False)
fig.tight_layout(rect=[0, 0.09, 1, 1]); fig.savefig("fig3_AB_curves.pdf"); plt.close(fig)
print("Fig 3: assertions passed (a(A)=lambda_h and a(B)<lambda_h for Thm 1; a(A)>=lambda_{h+1}, a(B)<lambda_{h+1} for Thm 3)")

# ---------------- numbers quoted in Remark 1 ----------------
gap = 1e9; arg = None
for h in range(2, 8):
    for N in range(4*h + 2, 4*h + 80):
        g = a_of(spider([h + 1, h, h, h] + [1]*(N - 2 - 4*h))) - lam(h + 1)
        if g < gap: gap, arg = g, (h, N)
print("Remark 1: min over 2<=h<=7, 4h+2<=N<=4h+79 of a(A)-lambda_{h+1} = %.3e at (h,N)=%s  (strictly positive: %s)" % (gap, arg, gap > 0))
