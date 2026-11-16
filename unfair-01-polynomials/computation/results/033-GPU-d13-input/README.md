# d13 necessary-prefix CPU preparation

This is a CPU-only preparation directory for the proposed numerical batch
`m=56`, `n=96`, first fractional index `d=13`, `starts=32768`,
`margin=.001`, and seed `2026111611`.  No command in this preparation started
MPS or a GPU worker.  The live d12 run was neither read for mutation nor
controlled by these scripts.

The complete raw ternary space has `3^12 = 531441` codes.  Integer arithmetic
finds 22,481 legal low prefixes and 12,291 survivors of the two declared
necessary filters.  Thus 32,768 rows cover the entire pool: 8,186 prefixes
receive three rows and 4,105 prefixes receive two.  No larger `starts` value
is proposed.

The copied worker is source SHA
`4cbe9ed5f04acebfaf6f420ce5a7c076b086352c317376cd14ef0520a8acff70`.
It differs from its frozen round032 source only by extending the explicit
audited `d` gate from 12 to 13 (both the public check and main assertion).
The independent generator and verifier do not import that worker.  The
verifier recomputes all 531,441 codes, checks every survivor, every one of the
32,768 cyclic assignment rows, all per-prefix counts, and then compares the
result to the worker's CPU-only audit.

Reproduce from this directory with:

```sh
/Users/zzh/.venvs/main/bin/python generate_d13_prefix_plan.py --out d13-prefix-plan.json --d 13 --m 56 --n 96 --requested-starts 32768 --margin .001 --seed 2026111611
/Users/zzh/.venvs/main/bin/python worker/gpu_factor_search_prefix.py --cpu-audit-only --m 56 --n 96 --fractional-index 13 --starts 32768 --steps 12000 --restart-every 2000 --margin .001 --seed 2026111611 --out worker-cpu-audit
/Users/zzh/.venvs/main/bin/python verify_d13_prefix_plan.py --input d13-prefix-plan.json --out independent-d13-plan-audit.json --worker worker/gpu_factor_search_prefix.py --prefix-module worker/fractional_prefixes.py --worker-cpu-audit worker-cpu-audit/cpu-prefix-audit.json
```

If the parent schedules this owned copy later, the exact bounded 4-hour
driver invocation is (not executed during this preparation):

```sh
/Users/zzh/.venvs/main/bin/python worker/run_prefix_batch.py --out live-m56n96d13 --m 56 --n 96 --d 13 --starts 32768 --steps 12000 --restart 2000 --margin .001 --seed 2026111611 --seconds 14400
```

The finite prefix conditions and the future floating-point trajectories are
necessary/numerical evidence only; they do not certify factorization or
exclusion.
