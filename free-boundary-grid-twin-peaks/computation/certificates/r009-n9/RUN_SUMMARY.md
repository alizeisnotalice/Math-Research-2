# n=9 extended-L completion kernel: corrected-source rerun

Status: a strict finite computation for the auxiliary completion kernel.  It
does not prove or disprove the original conditional two-peak distance problem.

## Source correction

The historical `boundary_kernel.cpp` hash-table builder used
`compare_exchange_weak`.  On a spurious failure, the expected value may remain
zero, after which the old code dereferenced `keys[e-1]`.  This rerun replaces
it with `compare_exchange_strong` and retains an explicit `e != 0` guard.
It also rejects a nonpositive thread count and asserts that stored indices fit
the uint32_t `index+1` slot representation and that the hash capacity does not
overflow.  The source was compiled without `-DNDEBUG`, because assertions
guard exact division and I/O operations.

The incomplete first attempt is deliberately retained as
`n9-r002-partial-090921/`: it created approximately 2.2 GB of segment files,
then its detached process disappeared without a diagnostic and before it
created any layer or result file.  It is not evidence for this calculation.
The authoritative rerun is `n9-r002/`, run in a live supervising process.

## Exact result

For the 9 by 9 extended-L prefix, the checked parameters are `S=17`,
`m=64`, boundary size `16`, `U=48`, and 14 ordered dependent pairs.  The
corrected rerun gives

\[
 h=2811307711161526712221155840260115808061564045287412733805364527470988492800,
\]
\[
 H=818140656621295511708762817109888719394876331671382960600810742737914675200,
\qquad T=H\cdot64\cdot63.
\]

Thus

\[
16H-4h=
1845019661294621338455581712717756278071765125592476434391513773922680832000>0,
\]

or `16 p(S) = 4.656285206336353... > 4`.  This reconfirms only the previously
identified failure of the auxiliary pointwise bound `p(S) <= 4/|boundary S|`.

## Independent checks

* Corrected C++ n=7 and n=8 reruns exactly match the old independent n=7
  Python-int result and n=8 Python certificate in `handoff-r9/`.
* `boundary_kernel_asan` was compiled with AddressSanitizer and UBSanitizer;
  the complete n=7 generate/calculate run has the exact expected output and
  no sanitizer diagnostic.  This is a small-instance memory-safety check, not
  a proof about every n=9 execution path.
* `independent_n9_layer_counts.py` does not read C++ mask layers.  A seven-row
  transfer DP gives the raw independent-set layer polynomial, while a separate
  reflection-quotient independent-set recurrence gives the diagonal-fixed
  polynomial.  It obtains raw total 877160008, fixed total 185378, and
  representative total 438672693.  These agree layer by layer with a fresh
  C++ audit of `n9-r002/layers`.

For the integer widths, every recurrence count is at most `64!`; a recurrence
numerator has at most 48 summands, so it is below `48*64! < 2^302`, within
U320.  The absolute alternating observable sum is at most
`877160008*14*64! < 2^330`, within U384.

## Reproduction and resources

Compile with:

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
c++ -O3 -std=c++20 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang \
  boundary_kernel.cpp -o boundary_kernel
```

Then run `boundary_kernel 9 12 n9-r002 generate`, followed by
`boundary_kernel 9 12 n9-r002 calculate`, and finally `audit`.  The actual
r002 generator took 10.15 seconds wall time.  The recurrence took 215.98
seconds wall, 518.46 user, and 124.64 system seconds; `/usr/bin/time -l`
reported 4,906,270,720 maximum resident bytes and zero swaps.  The system
peak-footprint field was 2,151,253,936 bytes; both raw tool fields are retained
in `n9-r002/calculate.stderr`.  The calculation process was absent after the
supervising process returned.

Files named `n9-result.json` and `n9-reflection-layer-audit-*.json` are n=9.
Historical bare names `handoff-r9/result.json` and
`handoff-r9/reflection-layer-audit.json` are n=8 files and must not be cited
as n=9 output.
