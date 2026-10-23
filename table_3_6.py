# Table 3.6: union of rectangles (Algorithm 3.2), omega* = 10, 20, start (0, 0)
import warnings
import numpy as np
from decimal import Decimal, ROUND_DOWN
import re
from f6_radius import N
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

A1 = np.array([[0.0, 0.1], [0.3, 0.0]])
for WY in (10, 20):
    G, NAME = 1.0, 'omega* = %d' % WY
    for ry in (1e-2, 1e-3):
        N[0] = 0 ; cols = column_union(A1, 0, 0, ry, G, WY)
        print(tr("%s r* = %g: %d columns, area %.4f, %d evaluations", (NAME, ry, len(cols), sum((c[1] - c[0]) * (c[3] - c[2]) for c in cols), N[0])))
        for c in cols: print(tr("   alpha in [%.6f, %.6f]  beta in [%.6f, %.6f]", c))
