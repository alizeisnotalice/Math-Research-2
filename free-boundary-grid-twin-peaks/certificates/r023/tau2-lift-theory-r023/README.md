# r023 tau=2 lift theory

`tau2-lift-theory.md` gives the exact first-merge conditional tail and second
moment of the stable-insertion gap \(L\), the all-size conditional bound, and
the resulting precise cross-fibre correlation gap.  It records sufficient
conditions only as unproved future targets.

Reproduce the finite exact n3 audit:

```sh
python3 /Volumes/data/CodexData/skills/python-math/scripts/run_math.py \
  --workdir /Volumes/data/CodexTemp/T2/a022-tau2-lift-theory-r023 \
  --log-root /Volumes/data/CodexTemp/T2/a022-tau2-lift-theory-r023/python-math \
  /Volumes/data/CodexTemp/T2/a022-tau2-lift-theory-r023/verify_tau2_lift_theory_n3.py
```

The script checks all 144 pairs \(1\le a,b\le12\) with exact fractions and
all \(9!\) words on the \(3\times3\) grid.  These checks are finite evidence;
the general statements in the note are proved algebraically.

The added `verify_fibre_moment_sharpness.py` checks (T7) on the exact parameter
grid \(1\le a,b\le100\), and evaluates the sharpness family
\(a=t,b=t^3\).  The analytic asymptotic in (T9a), rather than this finite
check, establishes that the constant \(4\) cannot be lowered for abstract
first-merge size parameters.
