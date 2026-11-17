# Actual public-package replay receipt

Replay was run in an isolated `_independent_replay/` directory on 2026-11-17.
The working source packages were not modified by the replay.

* Compiler: the recorded `short-gap-h45/reports/compiler-version.txt`; command
  was the C++20 command in `short-gap-h45/REPLAY.md`.
* A source smoke: compiled source run for `(n,v,prime,threads)=(3,4,0,1)`;
  resulting mathematical fields are saved in
  `short-gap-h45/reports/compiled-n3-v4-p0.json`. Its `h1`--`h5` residues
  are exactly `2688,1600,1024,640,384`, matching the saved checkpoint. The
  JSON files are not byte-identical only because `formula_status` changed from
  historical `pending independent proof` to final historical-status text; this
  is the documented I/O-only difference.
* A saved-checkpoint replay: `public-finite-crosscheck.json` reports 124
  checkpoint files, 31 points, a 252-bit CRT modulus, and `status: passed`.
  It actually read supplied r029 and r030 complete L references.
* B raw replay: gzip was decompressed in the isolated replay working directory;
  its SHA-256 was `9ec4cea2102d7ea1adb4b6757042e06e27ce375e0c53ce6c4a361664b9abe451`.
  The frozen f1593 auditor passed and every listed aggregate mathematical field
  matched the recorded aggregate within the verifier's `2e-13` tolerance.

The B result remains diagnostic only: `ESS=1.0000000003349663`,
`maxshare=0.9999999998325169`, and `Q_mRB=0.22509415843036934`. These do not
certify a probability, cone condition, mixing, SMC accuracy, or the original
mathematical question.
