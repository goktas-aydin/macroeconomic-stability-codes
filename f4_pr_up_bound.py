# Algorithm 2.3: upper bound u^*e of the extended omega*-Schur stability interval
import math
from numpy import linalg as LA
from f1_omega import omega

def pr_up_bound(A, B, ry, wy):
  # A: n x n numpy array with omega(A) <= wy; B: n x n numpy array (direction)
  # ry: r* > 0 (float); wy: omega* > 1 (float); returns the upper bound u^*e (float)
  n = A.shape[0] ; b = LA.norm(B, 2)
  i = 0 ; u= 0
  while True:
        i = i + 1
        w = omega(A) # Call omega function ; print("")
        if w is None: # Check if omega(A) returned None
            print("Omega could not be calculated for the current A matrix. Terminating up_bound.")
            return None # Indicate that up_bound also failed
        
        if w> wy:
            print("Not omega*-Schur stable (w > wy).")
            break

        a = LA.norm(A, 2) ;  r=(-a/b)+math.sqrt(a**2-(1/wy)+(1/w))/b
        if r < ry: break
        A = A + r * B ;  u = u+r ;
  return  float(u)
