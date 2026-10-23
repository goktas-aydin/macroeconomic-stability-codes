# Algorithm 3.1: single rectangle (edge policies 'alpha', 'sim', 'beta')
import numpy as np
from f6_radius import radius
from f7_edge_step import edge_step

E4 = {'R': (1, 2, 3, 'b', 1), 'L': (0, 2, 3, 'b', -1), 'T': (3, 0, 1, 'a', 1), 'B': (2, 0, 1, 'a', -1)}
def single_rectangle(A1, a0, b0, ry, g, wy, eps, policy):   # Algorithm 3.1, policy 'alpha', 'sim' or 'beta'
    # A1: 2 x 2 numpy array; (a0, b0): starting point (floats); ry: r* > 0; g, wy as in radius;
    # eps: first-stage threshold >= ry; policy: 'alpha', 'sim' or 'beta'
    # returns (alpha_lo, alpha_hi), (beta_lo, beta_hi)
    A0 = A1 + np.diag([a0, b0]) ; r0 = radius(A0, g, wy) ; R = [-r0, r0, -r0, r0]
    def sweep(e, tau):
        i, j, k, ax, d = E4[e] ; return edge_step(A0, R[i], R[j], R[k], ax, tau, g, wy)
    def cyclic(edges, tau):
        act = list(edges)
        while act:
            for e in list(act):
                s = sweep(e, tau)
                if s is None: act.remove(e)
                else: R[E4[e][0]] += E4[e][4] * s
    if policy == 'alpha': cyclic('RL', eps)
    if policy == 'beta': cyclic('TB', eps)
    while policy == 'sim':
        ss = [sweep(e, eps) for e in 'RTLB']
        if None in ss: break
        s = min(ss) ; R[:] = [R[0] - s, R[1] + s, R[2] - s, R[3] + s]
    def check():                                   # final check: one more sweep of each edge
        for e in order:
            s = sweep(e, ry)
            if s is not None: R[E4[e][0]] += E4[e][4] * s ; return True
        return False
    order = 'TRBL' if policy == 'beta' else 'RTLB'   # stage 2 order (Algorithm 3.1, step 5)
    cyclic(order, ry)
    while check(): cyclic(order, ry)
    return (a0 + R[0], a0 + R[1]), (b0 + R[2], b0 + R[3])
