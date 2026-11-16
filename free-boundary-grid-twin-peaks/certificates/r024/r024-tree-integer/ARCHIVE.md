# r024 independent integer tree DP

This package contains the direct rank-interleaving integer DP, separate from
the continuous-label message implementation. `run_verify.py` exhaustively
checks its small path/star cases and certifies the 30-vertex tree formed by
joining two height-three complete binary-tree roots. Its exact result is

\[
A_2/B=838907595/253657849>3=\Delta.
\]

Thus \(A_2\le\Delta B\) is false for all maximum-degree-three trees. It
does not refute that inequality for all degree-three *grids*, and does not
settle the original square-grid twin-peaks problem.
