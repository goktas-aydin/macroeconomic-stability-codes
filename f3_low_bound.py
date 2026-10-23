# Algorithm 2.2: lower bound l^e of the extended Schur stability interval
import math
from numpy import linalg as LA
from f1_omega import omega

def low_bound(A, B, ry, g):
    # A, B, ry, g as in up_bound; returns the lower bound l^e (float, <= 0)
    n = A.shape[0] ; b = LA.norm(B, 2)
    i = 0 ; l = 0
    while True:
        i = i + 1
        w = omega(A) # Call omega function
        if w is None: # Check if omega(A) returned None
            return None # Indicate that low_bound also failed
        
        a = LA.norm(A, 2)
        z = (-a / b) + math.sqrt(a**2 + (1 / w)) / b ;  r = z * g
        if r < ry: break
        A = A - r * B ;  l = l - r 
    return float(l)
