# r023 root activation scan

`root-activation.md` proves the exact low-root first-activation probability
under the uniform no-new-birth completion law and records the original finite
scan.  Its early candidate pointwise lower bounds are superseded by
`corner-redundant-root.md`: an insertion--deletion proof gives the sharp
all-state value \(\min p_{\rm act}=2/(n^2-k-1)\) for \(3\le k\le n^2-3\).
That theorem rules out the pointwise route only; it does not decide the
\(\tau=2\) weighted twin-peaks question.

Reproduce the exact integer n3/n4 scan:

```sh
python3 /Volumes/data/CodexData/skills/python-math/scripts/run_math.py \
  --workdir /Volumes/data/CodexTemp/T2/a022-root-activation-r023 \
  --log-root /Volumes/data/CodexTemp/T2/a022-root-activation-r023/python-math \
  /Volumes/data/CodexTemp/T2/a022-root-activation-r023/scan_root_activation.py
```

The enumerated universe is every nonempty connected vertex set \(C\) and
every vertex \(y\notin N[C]\), not merely states sampled from a two-peak
word.  The output stores exact completion counts and the individual neighbour
terms for each displayed minimizer.

`root_activation_exact.cpp` is an independent exact C++ scanner.  It validates
all n3/n4 extrema against the Python scan, then completed n5 in one process.
Its n5 \(2^{25}\)-entry `unsigned __int128` table is exactly 512 MiB; the
recorded run took 1.75 wall seconds, peaked at 538,591,232 resident bytes, and
reported zero swaps.  It uses no GPU.  Reproduce after checking current memory
pressure:

```sh
/usr/bin/clang++ -std=c++20 -O3 -DNDEBUG root_activation_exact.cpp -o root_activation_exact
/usr/bin/time -l ./root_activation_exact 5 n5-cpp.json
```

`corner_redundant_exact.cpp` fixes the corner/root subfamily independently.
The light-weight exact JSON cross-check (including every \(k\ge3\) for
\(n=3,4,5\)) is:

```sh
python3 /Volumes/data/CodexData/skills/python-math/scripts/run_math.py \
  --workdir /Volumes/data/CodexTemp/T2/a022-root-activation-r023 \
  --log-root /Volumes/data/CodexTemp/T2/a022-root-activation-r023/python-math \
  /Volumes/data/CodexTemp/T2/a022-root-activation-r023/verify_corner_redundant.py
```
