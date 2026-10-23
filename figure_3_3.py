# Figure 3.3: guaranteed regions, r* = 1e-3: union of rectangles (Algorithm 3.2) and single rectangle (Algorithm 3.1)
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from decimal import Decimal, ROUND_DOWN
import re
from f8_single_rectangle import single_rectangle
from f9_column_union import column_union

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

def draw(ax, cols, Ta, Tb, origin):
    for (a1, a2, b1, b2) in cols:
        if b2 - b1 < 1e-9: ax.plot([a1, a2], [b1, b1], color=BLUE_D, lw=3.0, solid_capstyle="butt")   # degenerate column
        else: ax.add_patch(Rectangle((a1, b1), a2 - a1, b2 - b1, facecolor=BLUE, edgecolor=BLUE_D, lw=0.8, alpha=0.9))
    ax.add_patch(Rectangle((Ta[0], Tb[0]), Ta[1] - Ta[0], Tb[1] - Tb[0], fill=False, edgecolor=ORANGE, lw=2))
    if origin: ax.plot([0], [0], "o", color=INK, ms=5, zorder=6)
    ax.plot([0.65], [0.7], "^", color=INK, ms=6, zorder=6)
    ax.set_aspect("equal") ; ax.set_xlabel(r"$\alpha$") ; ax.set_ylabel(r"$\beta$")
    ax.grid(color=GRID, lw=0.5) ; ax.set_axisbelow(True)

A1 = np.array([[0.0, 0.1], [0.3, 0.0]])
fig, axs = plt.subplots(2, 2, figsize=(6.3, 6.6))
for ax, g, wy, title in [(axs[0, 0], 0.95, None, r"(a) Schur mode, $\gamma=0.95$, start $(0,0)$"),
                         (axs[0, 1], 1.0, 10, r"(b) $\omega^*=10$, start $(0,0)$")]:
    cols = column_union(A1, 0, 0, 1e-3, g, wy) ; Ta, Tb = single_rectangle(A1, 0, 0, 1e-3, g, wy, 0.5, 'sim')
    draw(ax, cols, Ta, Tb, True) ; ax.set_title(title, loc="left") ; ax.set_xlim(-1.15, 1.15) ; ax.set_ylim(-1.15, 1.15)
cols = column_union(A1, 0.65, 0.7, 1e-3, 1.0, 10) ; Ta, Tb = single_rectangle(A1, 0.65, 0.7, 1e-3, 1.0, 10, 0.5, 'sim')
ax = axs[1, 0] ; draw(ax, cols, Ta, Tb, False)
ax.set_title(r"(c) $\omega^*=10$, start $(0.65,0.7)$", loc="left") ; ax.set_xlim(-1.15, 1.15) ; ax.set_ylim(-1.15, 1.15)
x0, x1, y0, y1 = 0.45, 0.85, 0.50, 0.90
ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor=INK, lw=0.8, ls=(0, (2, 2))))
ax = axs[1, 1] ; draw(ax, cols, Ta, Tb, False) ; ax.set_xlim(x0, x1) ; ax.set_ylim(y0, y1)
ax.set_title("(d) detail of (c): upward margins", loc="left")
ax.annotate("", xy=(Ta[1], 0.70), xytext=(0.65, 0.70), arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0))
ax.text((0.65 + Ta[1]) / 2, 0.69, tr(r"$%.3f$", Ta[1] - 0.65), ha="center", va="top", fontsize=8.5)
ax.annotate("", xy=(0.65, Tb[1]), xytext=(0.65, 0.70), arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0))
ax.text(0.645, (0.7 + Tb[1]) / 2, tr(r"$%.3f$", Tb[1] - 0.7), ha="right", va="center", fontsize=8.5)
print(tr("upward margins %.3f %.3f", (Ta[1] - 0.65, Tb[1] - 0.7)))
h = [Patch(facecolor=BLUE, edgecolor=BLUE_D, label="union of rectangles"),
     Line2D([0], [0], color=ORANGE, lw=2, label="single rectangle"),
     Line2D([0], [0], marker="o", color=INK, lw=0, label="starting point (0,0) in (a), (b)"),
     Line2D([0], [0], marker="^", color=INK, lw=0, label="reference point (0.65, 0.7); starting point in (c)")]
fig.legend(handles=h, loc="lower center", fontsize=7.5, frameon=False, ncol=2)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig("fig_3_3.pdf") ; fig.savefig("fig_3_3.png", dpi=130)
