# d15 MPS production manifest

- Launched: 2026-11-17 12:16 Asia/Shanghai; terminal session `35715`; OS PID
  `74343` (the PID is only an observation, not a completion claim).
- Command scope: numerical MPS exploration of d=15 prefixes, with `m=64`,
  `n=104`, `starts=49152`, `steps=12000`, `restart_every=2000`,
  `margin=0.001`, `seed=2026111701`, and `seconds=14400`.
- Live output: `d15-mps-production/`; config
  `d15-mps-production/config.json`; CPU input audit
  `d15-mps-production/cpu-prefix-audit.json`; periodic durable state
  `d15-mps-production/state.pt`.
- Production worker:
  `gpu_factor_search_d15.py`, SHA-256
  `5485fd6233a3dd62f1936dd20c23d2ab4effd9bec5187e5958cbcd8c91395de0`.
  It consumes the plan's saved survivor order, checks every consumed code with
  exact integer predicates, atomically writes JSON/checkpoints, and reports a
  normal completion or `partial_time_limit`.  It has no resume adapter yet;
  a partial checkpoint is not a claim of restartability.
- Plan: `d15-prefix-plan.json`, SHA-256
  `aa1d442efdb8efe1c1bd65a99f7ab8021222553008ecc915b13f9800323e5a7a`;
  4,782,969 raw codes, 94,971 valid low prefixes, 46,011 surviving the stated
  necessary filter, and 49,152 cyclic assignments.  The survivor-list SHA-256
  is `80f10e03d019aad3513165f940e45c311c50c3fb7dd981c3a4ab62cad0e009df`.
- Independent completed v1 plan replay is preserved at
  `verify_d15_prefix_plan_independent_v1_completed.py`, SHA-256
  `194f9fda058936becdb5d43e6e3e87cf808b5c587349d0f664bd814188493936`,
  with receipt `d15-prefix-independent-audit.json`.  It passed but its negative
  cases unnecessarily repeated enumeration; the later optimized source is not
  substituted for this completed audit.
- Full-size MPS smoke: `d15-mps-smoke/`, two updates in 4.868 seconds;
  maximum driver allocation 5,340,069,888 bytes.  Its `release.json` records
  `cleanup_called=true`, `process_will_exit=true`, and normal completion.
- Claim boundary: neither this finite prefix plan nor MPS residuals prove a
  global exclusion or furnish an exact factorization counterexample.
