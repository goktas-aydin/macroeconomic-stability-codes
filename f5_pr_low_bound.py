# Algorithm 2.4: lower bound l^*e of the extended omega*-Schur stability interval
import math
from numpy import linalg as LA
from f1_omega import omega

def pr_low_bound(A, B, ry, wy):
  # A, B, ry, wy as in pr_up_bound; returns the lower bound l^*e (float, <= 0)
  n = A.shape[0] ; b = LA.norm(B, 2)
  i = 0 ; l = 0
  while True:
        i = i + 1
        w = omega(A) # Call omega function
        if w is None: # Check if omega(A) returned None
            print("Omega could not be calculated for the current A matrix. Terminating low_bound.")
            return None # Indicate that low_bound also failed

        if w> wy:
            print("Not omega*-Schur stable (w > wy).")
            break
        a = LA.norm(A, 2) ; r=(-a/b)+math.sqrt(a**2-(1/wy)+(1/w))/b
        if r < ry: break
        A = A - r * B ;  l = l-r ;
  return  float(l)
