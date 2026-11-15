# R55 fixed-weight Cauchy screen: execution receipt

- MATLAB 26.2.0.3433840 (R2026b) Update 1; `-singleCompThread`; script-reported elapsed time 1.014685125 seconds.
- Script `cauchy_heavy_tail_fixed_screen.m`: SHA-256 `85b46ef12e307e72383876ad73cd29f583fb9616dcfba9f30ce40cac600272e1`.
- Input `input.json`: SHA-256 `2117f418133dae38994e361bb71d5dd6dff0a97d10e2a1f43027a074e6a8ad28`.
- Output `screen.json`: SHA-256 `83dd5605c1fb26fb012923fce00512f8cb9d17f8857edb76478e5a848966b6c3`.
- Working directory: `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-operations-terra/r055-v2`.
- Command: `/Applications/MATLAB_R2026b.app/bin/matlab -singleCompThread -batch "cauchy_heavy_tail_fixed_screen('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-operations-terra/r055-v2/input.json','/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-operations-terra/r055-v2')"`.
- Owned MATLAB PID 9880; wrapper observed exit status 0; subsequent `ps -p 9880 -o pid=,ppid=,command=` returned no process. No unrelated MATLAB process was stopped.
- An earlier run is excluded from reproducibility evidence: its source was changed during the run, so its executed version cannot be bound reliably.
- Floating screening only: finite grid and local minimization do not certify the global minimum, and numerical integration does not provide a rigorous enclosure.
