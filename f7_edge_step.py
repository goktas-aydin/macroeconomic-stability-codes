# edge step of Algorithm 3.1
import math
import numpy as np
from f6_radius import radius

def edge_step(A0, c, p, q, ax, tau, g, wy):   # edge {c}x[p,q] (ax='b') or [p,q]x{c} (ax='a')
    # A0: 2 x 2 numpy array at the starting point; c, p, q: floats, offsets from the starting point;
    # ax: 'b' (edge {c} x [p, q]) or 'a' (edge [p, q] x {c}); tau: threshold (float); g, wy as in radius
    # returns the advance s (float), or None if a step is smaller than tau
    t, s = p, math.inf
    while True:
        r = radius(A0 + np.diag([c, t] if ax == 'b' else [t, c]), g, wy)
        if r < tau: return None
        s = min(s, r)
        if t >= q: return s
        t = min(t + r, q)
