# omega(A): the Godunov-Bulgakov parameter ||H||, H from A H A^T - H + I = 0 (equivalent to (2.3) for the matrices of the paper, Section 2.2)
import numpy as np
from numpy import linalg as LA
from scipy import linalg

def omega(A):
  # A: real n x n matrix as a numpy array, e.g. np.array([[0.0, 0.1], [0.3, 0.7]])
  # returns omega(A) (float), or None if A is not Schur stable
  n = A.shape[0]; I = np.eye(n); w = None # Initialize w to None
  try:
    H = linalg.solve_discrete_lyapunov(A, I)
    if np.all(LA.eigvals(H) > 0): w = LA.norm(H, 2)
    else: print("The matrix A is not Schur stable.")
  except (linalg.LinAlgError, Exception): print("The matrix A is not Schur stable.")
  return w
