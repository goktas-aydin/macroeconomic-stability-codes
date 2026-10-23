# Example 3.16 and Table 3.7: start (0.65, 0.7), r* = 1e-3
import warnings
import numpy as np
from numpy import linalg as LA
from decimal import Decimal, ROUND_DOWN, ROUND_FLOOR, ROUND_CEILING
import re
from f1_omega import omega
from f6_radius import radius, N
from f8_single_rectangle import single_rectangle
from f9_column_union import column_union

def tr(fmt, *args, start=(0.0, 0.0)):   # as fmt % args, with the rounding of the paper:
    # %.nf: truncated (rounded towards zero); %.na, %.nb: alpha, beta rounded towards the starting point
    if len(args) == 1 and isinstance(args[0], tuple): args = args[0]
    it = iter(args)
    def rep(m):
        x = next(it)
        if m.group(2) is None: return m.group(0) % x
        k = m.group(2) ; x0 = 0.0 if k == 'f' else start[0 if k == 'a' else 1]
        mode = ROUND_DOWN if k == 'f' else (ROUND_FLOOR if x > x0 else ROUND_CEILING)
        return str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-int(m.group(1))), rounding=mode) + 0)
    return re.sub(r'%(?:\.(\d+)([fab])|[sdg])', rep, fmt)
warnings.filterwarnings('ignore')

A1 = np.array([[0.0, 0.1], [0.3, 0.0]])

At = A1 + np.diag([0.65, 0.7]) ; area = lambda cs: sum((c[1] - c[0]) * (c[3] - c[2]) for c in cs)
print(tr("||A~1|| = %.4f  omega(A~1) = %.4f", (LA.norm(At, 2), omega(At))))
for name, g, wy in [("Schur", 0.95, None), ("omega* = 20", 1.0, 20), ("omega* = 10", 1.0, 10)]:
    h = radius(At, g, wy) ; N[0] = 0 ; Ta, Tb = single_rectangle(A1, 0.65, 0.7, 1e-3, g, wy, 0.5, 'sim') ; n1 = N[0]
    N[0] = 0 ; cols = column_union(A1, 0.65, 0.7, 1e-3, g, wy)
    print(tr("%s: half-width %.4f  rectangle [%.4a, %.4a] x [%.4b, %.4b]  margins %.3f, %.3f  area %.3f  %d evaluations;  union %d columns, area %.3f, %d evaluations", (name, h, *Ta, *Tb, Ta[1] - 0.65, Tb[1] - 0.7, (Ta[1] - Ta[0]) * (Tb[1] - Tb[0]), n1, len(cols), area(cols), N[0]), start=(0.65, 0.7)))
for a, b in ((Ta[1], Tb[1]), (Ta[0], Tb[0])):          # corners as given (rounded towards (0.65, 0.7))
    ga, gb = float(tr("%.4a", a, start=(0.65, 0.7))), float(tr("%.4b", b, start=(0.65, 0.7)))
    print(tr("omega at the corner (%.4a, %.4b) = %.4f", (ga, gb, omega(A1 + np.diag([ga, gb]))), start=(0.65, 0.7)))
for c in cols: print(tr("   alpha in [%.4a, %.4a]  beta in [%.4b, %.4b]", c, start=(0.65, 0.7)))
cols = column_union(A1, 0.65, 0.7, 1e-2, 1.0, 10)
print(tr("omega* = 10, r* = 1e-2: union %d columns, area %.3f", (len(cols), area(cols))))
Ta, Tb = single_rectangle(A1, 0.65, 0.7, 1e-2, 1.0, 20, 0.5, 'sim')
print(tr("omega* = 20, r* = 1e-2: rectangle area %.3f", ((Ta[1] - Ta[0]) * (Tb[1] - Tb[0]))))
