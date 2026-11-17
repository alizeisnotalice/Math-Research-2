# Excluded initial background p0 attempt

The first `nohup` p0 process ended abnormally before it emitted a JSON result,
an exit-status receipt, or a `/usr/bin/time` receipt. It is not one of the
four canonical CRT residues. This package retains only its partial diagnostic
stderr as `raw/failed-nohup/n7-v24-p0.partial-stderr.txt` and records the
event for provenance.

There is no matching kernel, launchd, running-board, or jetsam diagnostic in
the inspected event window. Accordingly this archive does **not** call the
event OOM. The canonical p0 is the later foreground retry
`retry1-n7-v24-p0`, whose JSON, zero exit receipt, and resource receipt are in
`raw/canonical/` together with p1--p3.
