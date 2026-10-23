# Algorithm 3.2: union of rectangles
import numpy as np
from f6_radius import radius
from f7_edge_step import edge_step

def column_union(A1, a0, b0, ry, g, wy):                    # Algorithm 3.2
    # A1, a0, b0, ry, g, wy as in single_rectangle
    # returns the columns as a list of (alpha_lo, alpha_hi, beta_lo, beta_hi)
    A0 = A1 + np.diag([a0, b0]) ; xs = [0.0]
    for d in (1, -1):
        a = 0.0
        while (r := radius(A0 + np.diag([a, 0.0]), g, wy)) >= ry: a += d * r ; xs.append(a)
    xs = sorted(set(xs)) ; cols = []
    for lo, hi in zip(xs, xs[1:]):
        top = bot = 0.0
        while (s := edge_step(A0, top, lo, hi, 'a', ry, g, wy)) is not None: top += s
        while (s := edge_step(A0, bot, lo, hi, 'a', ry, g, wy)) is not None: bot -= s
        cols.append((a0 + lo, a0 + hi, b0 + bot, b0 + top))
    return cols
