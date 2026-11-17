# R20 local exact checks

All arithmetic is integer/GF(2), with no floating-point tolerance.

`mm_cases.py` builds293 standard-coordinate eight-variable MM pairs:288
four-changed-fiber cases,4 transposition-plus-full cases and1 two-full case.
The four-point domains use the two affine types {0,1,2,3} and {0,1,2,4};
all9 four-point derangements and all16 phase offsets are included. It checks
both phases by integer FWHT and each difference weight32. Python3.15.0;
initial measured elapsed0.039915041 s, exit0. Inputs bind this source hash.

`mm_cases_constructive_kernel_audit.py` independently reconstructs each phase
and constructs an explicit three-dimensional kernel using the proof's three
cases. It verifies the kernel avoids the deduplicated difference set. It
actually used12 workers and exited0. Its recorded script elapsed is 0.10319079200417036 s;
Python/platform and every worker PID are in the certificate. This is not a
scan of all97155 three-spaces, and does not address arbitrary bent pairs.

The root then checked every saved kernel independently by direct point/coset
collisions, without computing a difference mask. The persisted
`verify_local_certificates.py` was actually run and passed all293 certificates
and source/input bindings. It does not rerun the full FWHT audit.

`verify_signed_autocorrelation_obstacle_k2.py` constructs the real bent phase
f=x1*x2+x0*x3 and its legal flip on A={x3=0}. Both full Walsh transforms have
absolute value4. The signed partial transform has8 nonzero coefficients of
absolute value4. At delta=e1, signed autocorrelation is0 although8 support
pairs occur. This excludes replacing unsigned differences with signed
correlation zeros. It is not a counterexample to projectability: kernel<e3>
works. The root actually ran this final source; timing/version/hash are in
its JSON. No MATLAB was started; the finite local jobs exited.

Reproduce in a fresh data-volume temporary directory so archived outputs are
preserved. These scripts only require Python standard library. Run:

```sh
python3 mm_cases.py
python3 mm_cases_constructive_kernel_audit.py --input mm-cases.json --workers 12
python3 verify_signed_autocorrelation_obstacle_k2.py
python3 verify_local_certificates.py
```

The common-MM all-dimensional proof has two independent mathematical reviews.
Finite checks support the kernel constructions; they do not extend the
common-coordinate assumption to arbitrary bent pairs.
