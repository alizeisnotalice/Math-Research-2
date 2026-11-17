# Round 47 — completed local d15 window batch

This package freezes the completed (m=64,n=104,d=15) numerical window run
only.  It is not the subsequently started (m=80,n=120) run.

## Reproduction map

- `worker.py` is the exact production worker, SHA-256
  `d56d4b97dc691ccdab1c851d8e1ebb210c4941fc8e70b063d9c896c75ac77db2`.
- `frozen-base-worker.py` and `fractional_prefixes.py` are the byte-for-byte
  frozen dynamic dependencies, SHA-256 respectively
  `f0044bf5a0667b1825d1b58662243bc6c9b49682c809d045680ef5c53d348d9f`
  and `0b1a06770a61fe113031209986f676c7e1afb862d27c0c697dbebc0a9f78481d`.
- `prefix-plan.json` has SHA-256
  `aa1d442efdb8efe1c1bd65a99f7ab8021222553008ecc915b13f9800323e5a7a`;
  it contains 46,011 retained d15 codes.  The six production window offsets
  were 0, 8192, 16384, 24576, 32768 and 40960, so every retained code occurs
  in at least one window.
- `config.json` records the actual command parameters: MPS float32,
  `m=64,n=104,starts=8192,steps=6000,restart_every=1000,margin=.001`, and
  new seed `2026111702`.
- `cpu-prefix-audit.json`, `result.json`, and `release.json` are the actual
  output, numerical result, and exit receipt.  The process was PID 60800;
  `release.json` records normal completion, cleanup and a remaining MPS
  allocation of 31,687,424 bytes before process exit.

## Independent candidate status

The parent independently replayed all 16 saved candidate rows with its
Fraction validator and found no exact counterexample; the least exact decimal
residual was approximately `0.00326505908742547`.  This finite floating-point
batch neither excludes a factorization nor certifies a counterexample.


父独立十进制有理复核的可移植入口：

```sh
python3 verify_saved_samples.py --run . --verifier ../040-GPU-d14/verify_prefix_fraction.py --out /Volumes/data/CodexTemp/T2/round047-gpu-sample-check.json
```

这里只调用独立验证器的逐样本函数，不重新遍历4,782,969个三元前缀；完整前缀池的严格认证继承41轮已验输入。6×8192窗口的实际索引并集覆盖46,011个前缀，数值轨迹并不是逐前缀排除证明。
