# Algorithm 2.1: upper bound u^e of the extended Schur stability interval
import math
from numpy import linalg as LA
from f1_omega import omega

def up_bound(A, B, ry, g):
    # A: Schur stable n x n numpy array; B: n x n numpy array (direction), e.g. np.diag([1.0, 0.0])
    # ry: r* > 0 (float); g: gamma, 0 < g < 1 (float); returns the upper bound u^e (float)
    n = A.shape[0] ; b = LA.norm(B, 2)
    i = 0 ; u = 0
    while True:
        i = i + 1
        w = omega(A) # Call omega function
        if w is None: # Check if omega(A) returned None
            return None # Indicate that low_bound also failed
        
        a = LA.norm(A, 2)
        z = (-a / b) + math.sqrt(a**2 + (1 / w)) / b ;  r = z * g
        if r < ry: break
        A = A + r * B ;  u = u+r ;
    return  float(u)
