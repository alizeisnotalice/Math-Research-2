# n16 fixed-root production exit receipt

Command: `./new 16 119 65536 64 16 12 n16-root119-p65536-r64-l16-direct.csv 0xA0249000`.

The program exited with status 0.  The result CSV has its header plus 64 data
rows, one for each replicate 0 through 63.  A post-run process check found
none of the production shell, `/usr/bin/time`, or `new` process IDs
7943, 7944, and 7945 alive.

`/usr/bin/time -l` recorded 950.25 seconds wall time, 5301.00 seconds user
CPU time, 15.09 seconds system CPU time, maximum resident set size 421888000
bytes, peak footprint 422167272 bytes, and zero swaps.  These figures cover
this production process only; concurrent work belonging to other projects was
not stopped or attributed to it.

See `production-summary.json` for the strict input/output checks and finite
SMC diagnostic.  The CSV SHA-256 is recorded in `production-direct.sha256`.
