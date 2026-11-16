# 第 28 轮待整合：d=9 分数前缀 GPU 批次

本包封存完成的 `m=40,n=80,d=9` 分数前缀浮点 GPU 探索；它不是严格排除或精确反例。`queue.json` 记录 session 的 worker 正常完成（rc 0），而 `parent-exit-verification.json` 记录 session 34350、parent 13996 和 child 14000 已退出。

独立验证器不导入 worker 或 prefix module，精确枚举全部 `3^8=6,561` 个前缀：1,109 个低端合法，785 个通过两项已声明的必要过滤；SHA-256 分别为 `77f040bdd52553fd663a95256a3cc02a026f9b05edcfc20dec931bce2a9a9ed4` 与 `af61dace70e096bf15fdbd0514b52a7b505d778821a27520afc98d30ec1b833d`。每 cohort 的 32,768 starts 覆盖这 785 个池。

账本为五份 restart 快照（80 行）加 complete final cohort（16 行），共 96 个互异保存 JSON 样本；`partial=0`。每行以 `Fraction` 重卷积，复查 assigned/derived prefix、低卷积、d-cross、端点、非负及 margin。保存 residual 的最大差为 `4.73e-17`，最小保存的 Fraction residual 为 `0.0022413861006498337`；没有样本给出精确二元乘积。`release.json` 的退出前 allocator 值为 90,702,592 bytes，不能视为进程退出后的常驻内存。

保留 config/checkpoint/CPU prefix audit/release/queue/exit、5 份 restart、worker/driver/module、其引用的 frozen source，以及独立 verifier/receipt。排除可再生或冗余的 `state.pt`、`worker.log`、`result.json`。

复现：

```sh
python3 scripts/verify_prefix_fraction.py \
  --run-dir run \
  --worker scripts/gpu_factor_search_prefix.py \
  --prefix-module scripts/fractional_prefixes.py \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/prefix-d9-recheck.json
```
