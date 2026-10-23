# Figure 3.2: edge step of Algorithm 3.1 at the right edge of the initial square (Schur mode, gamma = 0.95)
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from decimal import Decimal, ROUND_DOWN
import re
from f6_radius import radius

def tr(fmt, *args):   # as fmt % args, but every %.nf value is truncated towards zero (as printed in the paper)
    if len(args) == 1 and isinstance(args[0], tuple): args = args[0]
    it = iter(args)
    def rep(m):
        x = next(it)
        if m.group(1) is None: return m.group(0) % x
        return str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-int(m.group(1))), rounding=ROUND_DOWN) + 0)
    return re.sub(r'%(?:\.(\d+)f|[sdg])', rep, fmt)
warnings.filterwarnings('ignore')

BLUE, BLUE_D, ORANGE, INK, GRID = "#86b6ef", "#2a78d6", "#eb6834", "#0b0b0b", "#efeeea"
# drawn at the printed size (text width 16 cm), so that all lettering is 7.5-8.5 pt (label R^(j) in Figure 3.2: 11 pt)
plt.rcParams.update({"font.family": "serif", "font.size": 8.5, "mathtext.fontset": "cm",
                     "axes.titlesize": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})

A1 = np.array([[0.0, 0.1], [0.3, 0.0]])
r0 = radius(A1, 0.95, None)                      # initial square of half-width r0 (Observation 1)
c, p, q = r0, -r0, r0                             # right edge {c} x [p, q]
nodes, rs, t = [], [], p
while True:
    r = radius(A1 + np.diag([c, t]), 0.95, None) ; nodes.append(t) ; rs.append(r)
    if t >= q: break
    t = min(t + r, q)
s = min(rs) ; print(tr("half-width %.6f  s = %.6f", (r0, s)))
fig, ax = plt.subplots(figsize=(4.55, 4.3))
ax.add_patch(Rectangle((-r0, -r0), 2 * r0, 2 * r0, facecolor=BLUE, edgecolor=BLUE_D, lw=1.0, alpha=0.55))
for i in range(len(nodes) - 1):
    ax.add_patch(Rectangle((c, nodes[i]), rs[i], nodes[i + 1] - nodes[i],
                           facecolor=("#f3f2ee" if i % 2 == 0 else "#e4e3df"), edgecolor=INK, lw=0.7))
k = 1
ax.add_patch(Rectangle((c - rs[k], nodes[k] - rs[k]), 2 * rs[k], 2 * rs[k], fill=False, edgecolor=INK, lw=0.8, ls=(0, (3, 2))))
ax.add_patch(Rectangle((c, p), s, q - p, facecolor=ORANGE, edgecolor=ORANGE, lw=1.4, alpha=0.45))
ax.plot([c] * len(nodes), nodes, "o", color=INK, ms=5, zorder=5)
for i, ti in enumerate(nodes):
    lab = r"$t_0=p$" if i == 0 else (r"$t_{m}=q$" if i == len(nodes) - 1 else r"$t_{%d}$" % i)
    ax.text(c - 0.045, ti, lab, ha="right", va="center", fontsize=8,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.85))
ax.annotate("", xy=(c + s, q + 0.07), xytext=(c, q + 0.07), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.6))
ax.text(c + s / 2, q + 0.11, r"$s=\min_i r_i$", ha="center", va="bottom", fontsize=8.5)
ym = nodes[k] + 0.5 * (nodes[k + 1] - nodes[k])
ax.annotate("", xy=(c + rs[k], ym), xytext=(c, ym), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.8))
ax.text(c + s + 0.5 * (rs[k] - s), ym + 0.02, r"$r_1$", ha="center", va="bottom", fontsize=8)
ax.text(-0.1, 0.0, r"$\mathcal{R}^{(j)}$", ha="center", va="center", fontsize=11)
ax.set_xlim(-0.75, 1.1) ; ax.set_ylim(-0.8, 0.95) ; ax.set_aspect("equal")
ax.set_xlabel(r"$\alpha$") ; ax.set_ylabel(r"$\beta$")
ax.grid(color=GRID, lw=0.5) ; ax.set_axisbelow(True)
h = [Patch(facecolor=BLUE, edgecolor=BLUE_D, alpha=0.55, label=r"guaranteed rectangle $\mathcal{R}^{(j)}$"),
     Line2D([0], [0], marker="o", color=INK, lw=0, label=r"nodes $(c,t_i)$ on the right edge"),
     Patch(facecolor="#e4e3df", edgecolor=INK, lw=0.7, label=r"rectangles $[c,c+r_i]\times[t_i,t_{i+1}]$"),
     Line2D([0], [0], color=INK, lw=0.8, ls=(0, (3, 2)), label=r"square of half-width $r_1$ (Theorem 2.2)"),
     Patch(facecolor=ORANGE, edgecolor=ORANGE, alpha=0.45, label=r"strip $[c,c+s]\times[p,q]$")]
fig.legend(handles=h, loc="lower center", fontsize=7.5, frameon=False, ncol=2)
fig.tight_layout(rect=(0, 0.16, 1, 1))
fig.savefig("fig_3_2.pdf") ; fig.savefig("fig_3_2.png", dpi=130)
