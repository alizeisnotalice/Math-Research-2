# First n=24 P17 trace attempt: retained failure record

This is a failure record, not a trace certification and not input to the n=24
production gate.

- Exact command arguments were `24 17 1 2 2 <first-raw> 9e270e0000000001 strat 16 trace` with the frozen binary.
- The first raw file has four JSONL lines (metadata, one replicate, one RB
  payload, run_end) and SHA-256
  `becde8798b957738c914ce82ed3020015c7991f5ffb74fcbe290ade84adbf308`.
- The expected sibling `<first-raw>.trace.json` was absent when the independent
  replay was invoked; that replay therefore raised `FileNotFoundError`.  The
  payload audit was not run for this first attempt.
- `/usr/bin/time -l` recorded 3.23 real seconds, 0 swaps and 0 received
  signals; its report SHA-256 is
  `806bd307d0f9240a7a63c16946df1dbe31e1dd5557df077a3e3e26bc88229456`.
  The captured stdout is empty (SHA-256
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`).
- The surrounding synchronous shell had `set -e` and advanced to the replay,
  so the wrapped `time` command returned zero.  No individual child PID was
  captured by that synchronous launch, hence no stronger per-PID exit receipt
  exists.  A later `ps` check found no live process referencing the first raw
  path, but that is only a later absence check.

Cause is **unestablished**.  In particular, this record makes no assertion
about environment flags, output-name behavior, or a program exception.  The
separate retry generated its own trace file and alone passed the rank-by-rank
and payload gates used for production.
