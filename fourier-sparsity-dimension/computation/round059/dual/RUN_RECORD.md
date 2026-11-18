# R59 positive weighted affine covers

For an ambient cube with dimension 3 or 4, each case is a family of distinct
nonempty affine flats. Repeated generated flats have their positive integer
weights added. The weighted indicators sum to the specified positive constant
at every point. Dividing by that constant gives a normalized cover. With
(m) distinct flats, (a) the maximum number of flats containing a single
point, (t) the binary rank of the sum of their normal spaces, and (r) the
real rank of their mean-centered indicators, the checks are (tle m-a)
and (tle r).

| Cases | Maximum activity greater than 1 | Equality (t=m-a) | Equality (t=r) |
| --- | ---: | ---: | ---: |
| Mac mini, 24,000,000 | 23,869,942 | 259,674 | 262,856 |
| MacBook, 30,000,000 | 29,837,814 | 326,021 | 329,840 |

C checks coverage by point evaluation and finds the binary rank by testing
annihilators of flat-direction bases. Metal checks coverage coefficient by
coefficient in the integer Walsh basis and finds binary rank directly from the
supports of the separate indicator spectra. Both use the integer matrix whose
entry for point (x
e 0) and flat (A_i) is (1_{A_i}(x)-1_{A_i}(0)).
Taking differences from one point is injective on the mean-zero subspace, so
this matrix has the claimed real rank. Entries are in ({-1,0,1}), and
each minor of order (kle15) has absolute determinant at most (k^{k/2} <
2147483647) by Hadamard. The rank over the prime field of that order is
therefore the exact real rank.

The indexed family consists of coordinate hyperplane controls, shuffled
unions of odd cycles, positive mixtures of affine partitions, and mixtures of
an odd-cycle cover and an affine partition. A cycle edge is a 2-point
affine flat. Within the disjoint union of odd cycles, giving each edge
weight (1/2) covers every vertex with total weight 1. Each odd component
rules out a disjoint edge partition. Controls 2 and 3 are explicit
instances. Controls 0 and 1 mix all coordinate hyperplane pairs with
unequal weights and reach equality in both bounds. The hash and
statistics match between CPU and Metal on each addressed range. A
separate Python calculation with exact rational elimination reproduces
the packed results of 16,384 initial case IDs and the four explicit
control certificates. Case IDs can generate the same cover.

The Mac mini and MacBook each used all 12 and 15 logical CPU worker slots in
parallel with Metal. The project CPU processes reached 784.4% and 902.7%
in the sampling window, where one occupied logical core represents 100%.
Their averages were 3.64 and 5.91 logical core equivalents respectively.
Both machines reached 100% whole-machine CPU and Metal device load. Whole-device
GPU load was already 99% and 100% before launch, and includes other work.
This run's GPU work is established by its own command-buffer receipts.
On both devices the CPU intervals intersected both sets of actual GPU
command intervals for an aggregate 0.980 seconds in 248 positive-length
intervals, as measured with host UTC clocks. Commands were issued to two
Metal queues per host. Summed GPU command durations therefore include
overlap. All four compute processes exited with status 0, with their PIDs
absent after exit.

## Reproduction

Both machines used Apple clang 21.0.0 and Apple Swift 6.4.2. The Mac mini
ran macOS 27.2 (26B123) on Apple M6; the MacBook ran macOS 27.3 (26C5035r)
on Apple M5 Pro. Python 3.9.6 supervised the local run. Python 3.15.0 ran the
independent integer and rational archive checks. Python 3.10 or newer is
required to regenerate the bank and run `verify_results.py`.

From this directory on an Apple Metal device, using available shared module
caches and appropriate temporary storage:

```sh
clang -O3 -std=c11 -pthread cover_cpu.c -o cover_cpu
swiftc -O CoverMetal.swift -o cover_metal
python3 -B make_bank.py bank-rebuilt.bin
cmp bank-rebuilt.bin bank.bin
python3 -B verify_results.py local macbook recreated-summary.json
cmp recreated-summary.json VERIFIED.json
shasum -a 256 -c SHA256SUMS
```

For the four-way run, copy the same sources and bank to both Apple Metal hosts.
Start `python3 -B run_host.py output 0 24000000` on the Mac mini and
`python3 -B run_host.py output 24000000 30000000` on the MacBook. Each
call requires a fresh output directory. After both write `ready.json`,
write a common available future Unix timestamp to each output directory's
`start_epoch.txt`. The archived records use 1794976129.464292. The
initial 16,384-case pilot independently used
`./cover_cpu bank.bin 0 16384 12 0 pilot_cpu.json pilot_cpu.csv` and
`./cover_metal bank.bin pilot_gpu.json 0 16384 8192`.
