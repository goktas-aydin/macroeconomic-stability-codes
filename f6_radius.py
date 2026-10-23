# step of Theorem 2.2 (gamma*delta(A)) or Theorem 2.3 (delta*(A)), (3.15)
import math
from numpy import linalg as LA
from f1_omega import omega

N = [0]                                   # number of evaluations of omega
def radius(A, g, wy):                     # g*delta(A) (Theorem 2.2, wy = None) or delta*(A) (Theorem 2.3)
    # A: 2 x 2 numpy array; g: gamma, 0 < g < 1 (Schur mode); wy: None (Schur mode) or omega* > 1
    # (then g is not used); returns the admissible step (float); counts the evaluations of omega in N[0]
    N[0] += 1 ; w = omega(A)
    if w is None or (wy and w > wy): return 0.0
    a = LA.norm(A, 2)
    if wy is None: return g * (-a + math.sqrt(a * a + 1 / w))
    return -a + math.sqrt(max(a * a + 1 / w - 1 / wy, a * a))
