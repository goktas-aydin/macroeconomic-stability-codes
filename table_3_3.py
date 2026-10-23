# Table 3.3: single rectangle (Algorithm 3.1), Schur mode, gamma = 0.95, start (0, 0)
import warnings
import numpy as np
from numpy import linalg as LA
from decimal import Decimal, ROUND_DOWN
import re
from f1_omega import omega
from f6_radius import radius, N
from f8_single_rectangle import single_rectangle

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
print(tr("||A1|| = %.4f  omega(A1) = %.4f  initial half-width r0 = %.6f (Schur, gamma = 0.95)", (LA.norm(A1, 2), omega(A1), radius(A1, 0.95, None))))
G, WY, NAME = 0.95, None, 'Schur'
for ry in (1e-2, 1e-3):
    for eps, pol in [(0.5, 'alpha'), (0.5, 'sim'), (0.5, 'beta'), (0.1, 'alpha'), (0.1, 'sim'), (0.1, 'beta')]:
        N[0] = 0 ; Ta, Tb = single_rectangle(A1, 0, 0, ry, G, WY, eps, pol)
        print(tr("%s r* = %g  eps = %g  %s  T_alpha = [%.6f, %.6f]  T_beta = [%.6f, %.6f]  area %.2f  %d evaluations", (NAME, ry, eps, pol, *Ta, *Tb, (Ta[1] - Ta[0]) * (Tb[1] - Tb[0]), N[0])))
