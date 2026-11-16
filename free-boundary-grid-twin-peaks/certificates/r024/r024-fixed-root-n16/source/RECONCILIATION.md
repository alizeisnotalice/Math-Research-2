# Fixed-root SMC audit scope (r024)

`fixed_root_n16.cpp` is a byte-for-byte copy of the candidate source at
`/Volumes/data/CodexTemp/T2/a022-fixed-root-n16-r024/fixed_root_n16.cpp`.
Its SHA-256 is `d8e7040f251c60ba7b660fcdba9e153b9dc21e351105fb474ad32314cc3456a4`.
The audit executable was compiled from that copy with
`clang++ -std=c++20 -O3 -pthread`; its SHA-256 is recorded in `audit.sha256`.

Against r022's `fixed_root_smc.cpp`, the candidate changes the capacity from
121 to 256 and replaces the `__uint128_t` masks by `std::bitset<256>`.  It
also changes the normalizer from the order count `h_root` to
`h_root/(N-1)!`: at rank `r` the factor is divided by `N-r+1`, and the product
of these denominators is `(N-1)!`.  The prefix proposal, systematic
resampling index construction, `mt19937_64` generator, unbiased bounded
integer routine, and endpoint-rejected move appear syntactically identical
apart from the bit-mask representation.

The old program has no terminal hash field.  Thus comparisons can prove the
same seed, normalizer conversion, minimum ESS, dead count, acceptance count,
and proposal count, but cannot prove identical terminal permutations.  The
new program independently checks every terminal particle with `peaks()==1`.
