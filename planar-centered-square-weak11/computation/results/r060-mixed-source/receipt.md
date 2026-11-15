# R60 mixed-time LP execution receipt

The frozen script SHA-256 is `db5edb77f26039c67aca344585d6476b3f9ba14e05950f2317a8b6199fbfac57`. It needs no input file. Its output `mixedtime-source-lp.json` has SHA-256 `fc8c1c457bc76cae8ad689c508e0784a3632dd3b12406b570ce8e105c270e88a`.

MATLAB 26.2.0.3433840 (R2026b) Update 1 ran with `-singleCompThread`; reported script time was 1.37763075 seconds. From this directory the command was `/Applications/MATLAB_R2026b.app/bin/matlab -singleCompThread -batch "cauchy_mixedtime_source_lp('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-operations-terra/r060-mixed-v2')"`. Owned PID 23923 exited 0; a subsequent `ps -p 23923 -o pid=,ppid=,command=` returned no process.

An earlier frozen script version failed during structure initialization, PID 22528 exit 1, and was verified exited. The successful LP solver reported ignoring 51 coefficients at or below 1e-9. The reported objective and full numerical evaluation are floating screens, not rigorous certificates.
