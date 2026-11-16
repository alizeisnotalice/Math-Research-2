# r024 fixed-root n=16 audit

The retained source extends the fixed-root kernel to 256 vertices. The audit
checks n=3/4 exact-root values, matching old n<=11 same-seed trajectories,
and independent n=16 terminal peak assertions. The production raw file is
the completed root 119=(7,7), P=65536, R=64, L=16 run; the receipt confirms
all 64 seeds, positive normalizers and terminal assertions. Run
`python3 -B verify_production_archive.py` for the retained source/raw/hash/
seed/ESS audit. Binaries, ASAN artifacts and redundant comparison CSVs are
excluded.

This is a finite one-peak SMC normalizer diagnostic, not an exact count.
