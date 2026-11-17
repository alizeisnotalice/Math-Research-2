# R20 reproduction index

See [local checks](RUN_RECORD_LOCAL.md), [MacBook CPU](RUN_RECORD_COMMON_MM_CPU.md),
and [two-host Metal](RUN_RECORD_METAL.md) for actual versions, seeds, commands,
timing boundaries, corrections and exit evidence. No finite computation is
claimed to settle general A28.

In a fresh data-volume temporary directory, copy this round's necessary files.
The C++ checker takes `smoke|full output.json`; compile with Apple clang21,
`-O3 -std=c++17 -pthread` and the shared data-volume compiler cache. The local
full run used the final agent-built binary whose SHA is in local-final-receipt.
Its actual program time was0.00525271s, exit0; the remote final run was
0.00695383s. Both use seed20261120. Standard-library shuffle implementation
is version-dependent; portable mathematical checking uses archived truth tables.

`python3 verify_local_certificates.py` checks the293 saved explicit kernels
by direct point collisions and verifies bindings. `python3 verify_cross_dimension.py`
checks both CPU receipts, equality of all120 case fields except timing, recovers
permutations/offsets from truth tables, and directly checks every kernel.
They are independent certificate checks, not additional full Walsh scans.
The final Metal source supports every N16..4096 through global staged buffers.

Final outputs and sources are hash-bound in manifest.json. The initial
zero-relative-offset sample was superseded by the actually rerun randomized
relative-offset sample; every nontrivial dimension/type stratum was independently
checked to include nonzero offsets. Local and remote CPU/GPU tasks all exited;
no sustained device saturation claim is made. No MATLAB session was started.
