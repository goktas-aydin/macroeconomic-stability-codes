# Figure 3.1: nodes of Algorithms 2.1-2.4 for system (3.4) and the step computed at each node
import math
import warnings
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from decimal import Decimal, ROUND_DOWN
import re
from f1_omega import omega

def tr(fmt, *args):   # as fmt % args, but every %.nf value is truncated towards zero (as printed in the paper)
    if len(args) == 1 and isinstance(args[0], tuple): args = args[0]
    it = iter(args)
    def rep(m):
        x = next(it)
        if m.group(1) is None: return m.group(0) % x
        return str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-int(m.group(1))), rounding=ROUND_DOWN) + 0)
    return re.sub(r'%(?:\.(\d+)f|[sdg])', rep, fmt)
warnings.filterwarnings('ignore')

def nodes(A, B, d, g, wy):     # nodes and steps of Algorithms 2.1-2.4 without the stopping test
    # A, B: n x n numpy arrays; d: +1 (upper bound) or -1 (lower bound); g: gamma (Schur mode);
    # wy: None (Schur mode) or omega*; continues until the step falls below 1e-5
    # returns the nodes xs and the steps st computed at them (st[k] at xs[k])
    b = LA.norm(B, 2) ; xs, st, x = [0.0], [], 0.0
    while True:
        w = omega(A)
        if w is None or (wy is not None and w > wy): st.append(0.0) ; break
        a = LA.norm(A, 2)
        r = g * ((-a / b) + math.sqrt(a ** 2 + 1 / w) / b) if wy is None else (-a / b) + math.sqrt(a ** 2 - 1 / wy + 1 / w) / b
        st.append(r)
        if r < 1e-5: break
        A = A + d * r * B ; x = x + d * r ; xs.append(x)
    return xs, st

def endpoint(xs, st, ry):       # endpoint of the stopping rule: stop at the first step smaller than r*
    k = 0
    while k < len(st) and st[k] >= ry: k += 1
    return xs[k] if k < len(xs) else xs[-1]

BLUE, BLUE_D, ORANGE, INK, GRID = "#86b6ef", "#2a78d6", "#eb6834", "#0b0b0b", "#efeeea"
# drawn at the printed size (text width 16 cm), so that all lettering is 7.5-8.5 pt (label R^(j) in Figure 3.2: 11 pt)
plt.rcParams.update({"font.family": "serif", "font.size": 8.5, "mathtext.fontset": "cm",
                     "axes.titlesize": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})

Ab, E11 = np.array([[0.0, 0.1], [0.3, 0.7]]), np.diag([1.0, 0.0])
fig, axs = plt.subplots(1, 2, figsize=(6.3, 3.25), sharey=True)
for ax, wstar, title in zip(axs, [None, 10], [r"(a) Schur mode, $\gamma=0.95$", r"(b) $\omega^*$ mode, $\omega^*=10$"]):
    up = nodes(Ab, E11, +1, 0.95, wstar) ; dn = nodes(Ab, E11, -1, 0.95, wstar)
    for xs, st in (up, dn):
        ax.plot(xs[:len(st)], st, "-", color=BLUE_D, lw=1.0, zorder=2)
        ax.plot(xs[:len(st)], st, "o", color=BLUE_D, ms=2.6, zorder=3)
    for ry, ls in zip([0.1, 0.01, 0.001], ["-", "--", ":"]):
        lo, hi = -abs(endpoint(*dn, ry)), endpoint(*up, ry)
        ax.axhline(ry, color=ORANGE, lw=1.2, ls=ls, zorder=1)
        ax.text(-0.42, ry * 1.12, r"$r^*=%g$" % ry, fontsize=7.5, va="bottom", ha="center")
        ax.plot([lo, hi], [ry, ry], "|", color=INK, ms=8, mew=1.3, zorder=4)
        print(tr("%s  r* = %g  endpoints %.6f, %.6f", ("Schur" if wstar is None else "omega* = 10", ry, lo, hi)))
    hi = endpoint(*up, 0.1)
    name = r"u^{e}(0.1,0.95)" if wstar is None else r"u^{*e}(0.1,10)"
    src = "Corollary 3.2" if wstar is None else "Table 3.1"
    txt = tr(r"$%s=%.5f$" if wstar is None else r"$%s=%.6f$", (name, hi))
    ax.annotate(txt + "\n(" + src + ")", xy=(hi, 0.1), xytext=(hi - 0.95, 0.0035), fontsize=7.5,
                ha="left", arrowprops=dict(arrowstyle="->", color=INK, lw=0.8))
    ax.set_yscale("log") ; ax.set_xlim(-1.06, 1.0) ; ax.set_ylim(5e-5, 0.7)
    ax.set_xlabel(r"node $u_k$, $l_k$ (value of $\alpha$)") ; ax.set_title(title, loc="left")
    ax.grid(color=GRID, lw=0.5, which="both") ; ax.set_axisbelow(True)
axs[0].set_ylabel("step computed at the node")
h = [Line2D([0], [0], marker="o", color=BLUE_D, lw=1.0, ms=3, label="node and the step computed there"),
     Line2D([0], [0], color=ORANGE, lw=1.2, label=r"threshold $r^*$"),
     Line2D([0], [0], marker="|", color=INK, lw=0, ms=8, mew=1.3, label="endpoint")]
fig.legend(handles=h, loc="lower center", fontsize=7.5, frameon=False, ncol=3)
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig("fig_3_1.pdf") ; fig.savefig("fig_3_1.png", dpi=130)
