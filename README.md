# macroeconomic-stability-codes

Python codes for R. N. Göktaş and K. Aydın, *Sensitivity and Guaranteed Schur Stability in a
Two-Sector Macroeconomic Model: A Computational Godunov–Bulgakov Approach.*

## Usage

In this folder (Python 3.8 or later):

```
pip install -r requirements.txt
python table_3_3.py
```

## Files

- `f1_omega.py` – `f9_column_union.py`: ω(A) and Algorithms 2.1–2.4, 3.1, 3.2
  (`ry` = r\*, `g` = γ, `wy` = ω\*; `wy = None` in the Schur mode)
- `table_3_3.py` – `table_3_6.py`, `example_3_16.py`: Tables 3.3–3.7 and Example 3.16
- `figure_3_1.py` – `figure_3_3.py`: Figures 3.1–3.3

## Other results

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
| Table 2.1 | `omega(np.array([[-z, 10**(k-1) - 1], [0, z]]))`, z = 0, 0.1, k = 1, …, 5 |
| Theorem 3.1 | T = (−u, u), `u = -a + math.sqrt(a**2 + 1/omega(A1))`, `a = LA.norm(A1, 2)` |
| Table 3.1, Part II | `pr_low_bound(A1, E11, 0.1, wy)`, `pr_up_bound(A1, E11, 0.1, wy)`, wy = 20, 30, 40 |
| Example after Theorem 3.6 | `B = np.array([[0.65, 0.1], [0.3, 0.0]])`; `low_bound(B, E22, 0.1, 0.95)`, `up_bound(B, E22, 0.1, 0.95)`, `pr_low_bound(B, E22, 0.1, 10)`, `pr_up_bound(B, E22, 0.1, 10)` |
| Remark 3.8(i) | `low_bound(A1, E11, 0.1, g)`, g = 0.855, 0.857; nodes: loop below |
| Table 3.2 | `up_bound(A1, E11, 1e-5, 1 - 10**-k)`, `pr_up_bound(A1, E11, 1e-5, 10**(2*k - 1))`, k = 1, …, 5 |

```python
for g in (0.855, 0.857):
    A, l = A1.copy(), 0.0
    while True:
        a = LA.norm(A, 2) ; r = g * (-a + math.sqrt(a**2 + 1 / omega(A)))
        print(g, l, r)
        if r < 0.1: break
        A, l = A - r * E11, l - r
```

## Notes

The values are given rounded towards the starting point, as in the paper (truncated for the starting
points 0 and (0, 0)). The paper uses Python 3.9.12, NumPy 1.21.5, SciPy 1.7.3 and Matplotlib 3.5.1;
current versions give the same printed values (the figures differ only slightly in lettering).

## License

BSD 3-Clause (`LICENSE`).
