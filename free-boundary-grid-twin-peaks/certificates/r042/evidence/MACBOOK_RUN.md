# Remote exact curved-cut run receipt (r073)

Remote task: `/Users/zhengzhihao/CodexRemoteMath/a022-curved-cut-second-moment-r073`.

## Frozen inputs and gate

- Curved-cut source SHA-256: `4258a7947abe18bbadb28f9d48d7ac87403ea8e3b5c3508a6d67d37a632786c5`
- Remote Apple-clang binary SHA-256: `618ef2511e00d8f69cb7e6a90e20dc24b3c038d094b78d465a8e46af7bb76f13`
- Orbit manifest SHA-256: `52969f98715c342d31ef6000c395316ee643b8d6d4b226fc38dc310c66d4a077`
- New/old exact-DP gates agreed on n=4 `2222`, n=5 `22222`, n=5
  `33333`, and their intersection. The corresponding counts were
  `41328967680`, `301412212232374763520`,
  `301412212232374763520`, and `14001081747604852400`.

The first detached launch did not start a worker because the remote reports
directory was absent before shell redirection. `gate/launch_attempt1_no_start.json`
records that zero-start event. The production launch then used fresh `raw/` and
`receipts/` directories.

## Production and local independent audit

The remote controller ran 152 symmetry-orbit jobs, at a configured cap of eight
workers. It recorded a zero exit code, validated count, output SHA-256, resource
SHA-256, and wait/reap receipt for every job. Its `final.json` records
`completed=152`, no failed or unlaunched jobs, and
`all_children_waited_reaped=true`.

The local independent audit read all 152 raw outputs, receipts, and resource
files; checked frozen hashes, output/receipt hash links, curve identities,
nonnegative U128 counts, and the diagonal/off-diagonal multiplicities
`32`/`496`. It obtained

\[
S_1=5596074816724757668240,
\qquad S_2=46414052741526275088808.
\]

With the supplied exact \(A_2(5)=7063261481362958280000\), the finite
second-moment lower bound is

\[
\frac{S_1^2}{A_2(5)S_2}
=\frac{4341733217943243972716111325054442367}
{45451778910155744334389011477607669550}
\approx0.09552394476188757514.
\]

The greatest recorded per-job RSS was 538,525,696 bytes. At 2026-11-17T14:11:51Z,
the remote named-task process count was zero and both the controller (99376) and
its caffeinate process (99377) were absent.

This is a strict finite \(n=5\) computation for the stated event family. It is
not a positive-probability subsequence result and does not settle the original
asymptotic twin-peaks question.
