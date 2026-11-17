# Rejected debug output

`results.json` produced before 2026-11-17 10:45 Asia/Shanghai is not a
mathematical search result and must not be cited.  Its extension predicate
incorrectly required every XOR of selected affine-flat indicators to itself be
one affine-flat indicator.  A legal Fourier sign mask can be a union/XOR of
flats, so this rejects the elementary eight-slice MM family `v=a`.

The corrected program uses a cached exact integer Walsh legality oracle for
every arbitrary XOR mask.  The retained file is debugging evidence only; no
negative conclusion, upper bound, or coverage claim follows from it.
