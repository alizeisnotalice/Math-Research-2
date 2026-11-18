# MacBook scan of the first three periodic output harmonics

For a (2pi)-periodic variance profile (a), put
[
 C_n=
rac1{2pi}int_0^{2pi}
       e^{in\theta-n^2 a(\theta)/2},d\theta,qquad
 Z_n=
rac1{2pi}int_0^{2pi}e^{-n^2 a(\theta)/2},d\theta.
]
The score is (max_{1le nle3}|C_n|/Z_n). It has the same
zeros as the three raw harmonics and is invariant under adding
a constant to (a). This prevents a large positive variance offset
from making a profile appear to cancel. The recorded coefficients
give a nonnegative variance profile after adding the sum of
individual mode amplitudes. The integrals use radians.

The reference part of the profile is drawn from five families,
each stratified through eight nominal coefficient norms
(1/64,1/32,1/16,1/8,1/4,1/2,3/4,1).

- **Low only** uses independently drawn sine and cosine modes
  (1,2,3), scaled to the layer norm.
- **High only** uses (4,5,7) at the layer norm.
- **Second-order corrected** starts with the same high-mode
  ensemble and adds low modes centered at their second-order
  cancellation formula. The randomized coordinate deviation
  is at most (0.75\alpha^3).
- **Subperiod control** uses modes (4,8), for which the first
  three harmonics vanish exactly by shorter-period symmetry.
- **Mixed** starts with modes (4,5,7) and adds uncorrected low
  coordinates with magnitude up to (0.3\alpha).

Here (alpha) is the norm assigned before any low-mode
augmentation. The high-only family does not impose a lower
bound on each active mode amplitude. Some low scores therefore
approach the single-mode symmetry.

## What the first three zero conditions allow

At a fixed positive constant (a=c), regard the sine and cosine
coefficients of modes (1,2,3) as six real unknowns. The derivative
of the six real harmonic equations is block diagonal, with
(n)-th block (-n^2 e^{-n^2c/2}I_2/4). The implicit function
theorem gives exact simultaneous zeros of these first three
harmonics for all sufficiently small high-mode coefficients.
Choosing nonzero high coefficients in modes 4 and 5 gives genuine
fundamental period (2pi), and positivity persists near the
constant. The fourth harmonic is nonzero for sufficiently small
such profiles with a nonzero first-order mode-4 coefficient.
The proof only uses this finite system.

Expanding about the constant yields the initial correction
(a_n=(n^2/4)(h^2)_n+O(|h|^3)) for (1le nle3), where
(h) is the higher-mode component and the coefficients on both
sides use the same cosine/sine convention. This formula seeds
one of the numerical families.

## Results and independent verification

The sweep used 512 periodic trapezoidal nodes and tested
**1,895,301,120** profiles across 40 strata. Each stratum had
47,382,528 profiles. All observed quantities were finite.

| Family | Best sampled score at (alpha=1/64) | Count below (10^{-5}) | Best sampled score at (alpha=1) |
| --- | ---: | ---: | ---: |
| Low only | 0.0037565245 | 0 | 0.16332328 |
| High only | (9.619429\times10^{-7}) | 10,087 | 0.0019155336 |
| Second-order corrected | (1.5175753\times10^{-7}) | 44,997,064 | 0.029070327 |
| Subperiod control | (1.6184458\times10^{-9}) | 47,382,528 | (2.7957165\times10^{-9}) |
| Mixed | 0.0003093244 | 0 | 0.011053908 |

The residuals in the subperiod control are Float32 errors on
structurally zero quantities. For each of the 40 minima, Float64
checks at 2048 and 4096 nodes agreed to better than (10^{-12}).
The largest absolute discrepancy from the recorded Metal
coordinates was (2.36237\times10^{-7}).

With the higher coefficients of each of the eight selected
second-order-corrected candidates fixed, a Float64 fit of six
low coordinates attained maximum residual below
(9.32\times10^{-16}) for every layer, including (alpha=1).
These fits provide numerical evidence at the selected profiles.

For a separate fixed profile with high coefficients
(c_4=0.05,c_5=0.04,s_7=0.02), both Macs solved the first six
coordinates with 86-digit arithmetic at 384 nodes and checked
them at 768 nodes. The largest absolute normalized coordinate
among the first three harmonics was about
(2.88034\times10^{-87}). The fourth normalized magnitude
was (0.1951239139584). This computation has finite precision
and accompanies the local exact implicit-function argument.

## Reproduction and shutdown

The MacBook used Apple M5 Pro with macOS 27.3 build 26C5035r and
Swift 6.4.2. The main Metal PID 87646 ran 03:19:28–03:21:55
UTC on 2026-11-18 and exited zero. Its command buffers accumulated
**100.018899 seconds**. The 227 GPU-device readings were 99–100%,
with mean 99.996%, and the initial device reading was already
100%. A post-run process query found neither the main nor the
smoke-test PID.

From the authenticated MacBook directory
`~/CodexRemoteMath/cc63-metal-nov18`:

```sh
export TMPDIR="$HOME/CodexRemoteMath/cache/cc63-temp"
xcrun swiftc -O -module-cache-path "$HOME/CodexRemoteMath/cache/swift" \
  -Xcc -fmodules-cache-path="$HOME/CodexRemoteMath/cache/clang" scan.swift -o harmonics_scan
sh run.sh smoke 256 0 1024 918
sh run.sh sync 512 100 65536 3719
"$HOME/.venvs/math/bin/python" check_cpu.py sync.json sync.cpu.json
```

The time-limited native run stops after its first full stratum
cycle over the target GPU command duration. The seed fixes the
candidate stream and a rerun may complete a different number
of cycles. The production run completed 723 cycles.
The MacBook CPU verification ran Python 3.12.13, SciPy 1.18.1,
and mpmath 1.3.0. The local verification produced identical
numerical output. Selected sources, statistics, and telemetry
match the MacBook hashes in `remote.sha256`.
