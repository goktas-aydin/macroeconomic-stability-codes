# macroeconomic-stability-codes

Codes for verifying the results of R. N. Göktaş and K. Aydın, *Sensitivity and Guaranteed Schur Stability
in a Two-Sector Macroeconomic Model: A Computational Godunov–Bulgakov Approach.*

## Functions

1. `f1_omega.py` – `omega(A)`, the parameter ω(A) (Section 2.2)
2. `f2_up_bound.py` – `up_bound(A, B, ry, g)`, Algorithm 2.1
3. `f3_low_bound.py` – `low_bound(A, B, ry, g)`, Algorithm 2.2
4. `f4_pr_up_bound.py` – `pr_up_bound(A, B, ry, wy)`, Algorithm 2.3
5. `f5_pr_low_bound.py` – `pr_low_bound(A, B, ry, wy)`, Algorithm 2.4
6. `f6_radius.py` – `radius(A, g, wy)`, step of Theorem 2.2 / 2.3, (3.15)
7. `f7_edge_step.py` – `edge_step(A0, c, p, q, ax, tau, g, wy)`, edge step of Algorithm 3.1
8. `f8_single_rectangle.py` – `single_rectangle(A1, a0, b0, ry, g, wy, eps, policy)`, Algorithm 3.1
   (edge policies 'alpha', 'sim', 'beta')
9. `f9_column_union.py` – `column_union(A1, a0, b0, ry, g, wy)`, Algorithm 3.2

Here `ry` is r\*, `g` is γ and `wy` is ω\* (`wy = None` in the Schur mode).

## Scripts

Run the scripts in this folder (they import the functions above), e.g. `python table_3_3.py`.

- `table_3_3.py` – Table 3.3 and Observations 1, 2, 4 (Schur mode; all edge policies)
- `table_3_4.py` – Table 3.4 and Observations 1, 2, 4 (ω\* = 10, 20; all edge policies)
- `table_3_5.py` – Table 3.5
- `table_3_6.py` – Table 3.6
- `example_3_16.py` – Example 3.16 and Table 3.7
- `figure_3_1.py` – Figure 3.1 (`fig_3_1.pdf`); also prints the endpoints of the examples after
  Corollaries 3.2 and 3.4 and of Part I of Table 3.1
- `figure_3_2.py` – Figure 3.2 (`fig_3_2.pdf`); also prints the half-width and s of its caption
- `figure_3_3.py` – Figure 3.3 (`fig_3_3.pdf`); also prints the upward margins of panel (d)

## Verifying the remaining results

The remaining values follow from single calls of the functions. These calls start at 0, so the paper gives their values truncated (rounded towards zero) to the digits shown.
In Python, in this folder:

```python
import math, numpy as np
from numpy import linalg as LA
from f1_omega import omega
from f2_up_bound import up_bound
from f3_low_bound import low_bound
from f4_pr_up_bound import pr_up_bound
from f5_pr_low_bound import pr_low_bound
A1 = np.array([[0.0, 0.1], [0.3, 0.7]]) ; E11 = np.diag([1.0, 0.0]) ; E22 = np.diag([0.0, 1.0])
```

| Result | Call |
|---|---|
| Table 2.1 | `omega(np.array([[-z, 10**(k-1) - 1], [0, z]]))` for z = 0, 0.1 and k = 1, …, 5 (for k = 5 SciPy may warn that the Lyapunov system is ill-conditioned, since ω is very large) |
| Theorem 3.1 | `a = LA.norm(A1, 2)`, `w = omega(A1)`; T = (−u, u) with `u = -a + math.sqrt(a**2 + 1/w)` |
| Table 3.1, Part II | `pr_low_bound(A1, E11, 0.1, wy)`, `pr_up_bound(A1, E11, 0.1, wy)` for wy = 20, 30, 40 |
| Example after Theorem 3.6 | with `B = np.array([[0.65, 0.1], [0.3, 0.0]])`: `low_bound(B, E22, 0.1, 0.95)`, `up_bound(B, E22, 0.1, 0.95)`, `pr_low_bound(B, E22, 0.1, 10)`, `pr_up_bound(B, E22, 0.1, 10)` |
| Remark 3.8(i) | `low_bound(A1, E11, 0.1, g)` for g = 0.855, 0.857; nodes and steps: see the loop below |
| Table 3.2 | `up_bound(A1, E11, 1e-5, 1 - 10**-k)`, `pr_up_bound(A1, E11, 1e-5, 10**(2*k - 1))` for k = 1, …, 5 |

```python
for g in (0.855, 0.857):          # nodes and steps of Algorithm 2.2 in Remark 3.8(i)
    A, l = A1.copy(), 0.0
    while True:
        a = LA.norm(A, 2) ; r = g * (-a + math.sqrt(a**2 + 1 / omega(A)))
        print(g, l, r)
        if r < 0.1: break
        A, l = A - r * E11, l - r
```

## Versions

All results of the paper were computed with Python 3.9.12, NumPy 1.21.5, SciPy 1.7.3 and Matplotlib 3.5.1
(Anaconda distribution); with these versions the figure scripts reproduce Figures 3.1–3.3 exactly. The
scripts run unchanged and print the same values with Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1 and
Matplotlib 3.11.2 (the figures then differ only slightly in lettering and layout). Python 3.8 or later is
required (`f9_column_union.py` uses `:=`).

## Rounding

As in the paper, the scripts give the computed values rounded towards the starting point of the
computation (function `tr`): values of α and β towards the corresponding coordinate of the starting point,
all other values towards zero. For the starting points 0 and (0, 0) this is truncation; only
`example_3_16.py` and panels (c), (d) of `figure_3_3.py` start elsewhere, at (0.65, 0.7) (`figure_3_3.py`
prints only the margins, which are truncated).

## License

BSD 3-Clause (`LICENSE`).
