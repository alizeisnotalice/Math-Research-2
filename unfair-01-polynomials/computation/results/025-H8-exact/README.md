# H8(x²)H8(x³) 的完整精确证书

固定 H8=1+x+x²+x⁷+x⁸。一般均匀和猜想仍开放。此包严格排除固定目标的所有非标准归一化非负分解，恒定分解另行处理。

从仓库根运行（Python 标准库即可）：

```sh
cd unfair-01-polynomials/computation/results/025-H8-exact
python3 verify_H8_exact.py --input H8-a2b3-inputs.json --balls H8-a2b3-root-balls.json --certificate H8-a2b3-certificate.json --galois H8-mother-certificate.json --record verification-reproduced.json
```

验证器独立生成全部 1482 种布局，不读搜索器 options 缓存。它重新验证模素数分解、8 个母根及 40 个目标根的有理 Rouché 圆盘、圆盘互异与共轭对应；用 Fraction 区间卷积复算每个被指定的负系数。1481 个布局的系数区间上界严格小于零，最小认证负余量为 102346591/137438953472。余下 E=6/code=0 精确为 H8(x²)/H8(x³)。因子的首一性来自原模型端点与系数上界，没有手动把数值端点改为 1。

根布局覆盖的解析依赖见 ../../../proofs 的第25轮 H8 稿及第18轮根幂积平衡引理。最终链接以本轮 proof 文件为准。算术证书的检查不使用浮点误差阈值。输入中的历史绝对路径只记来源，不是验证所需的运行时依赖。

搜索和验证分离：巨大完整区间中间件可重生成，公开包仅存外扩后的二进有理区间。compact certificate 的原始来源哈希用于溯源；每个区间仍由本包从根圆盘重新计算验证，不信任搜索成功标记。

生成链可重现如下（中间 raw 文件请写到 CodexTemp）：

```sh
python3 generate_h8_balls.py --input H8-a2b3-inputs.json --output /Volumes/data/CodexTemp/T2/unfair-01-polynomials/reproduced-balls.json
python3 search_h8_negative_intervals.py --input H8-a2b3-inputs.json --options H8-a2b3-options.json --balls H8-a2b3-root-balls.json --output /Volumes/data/CodexTemp/T2/unfair-01-polynomials/reproduced-raw.json
python3 compact_h8_intervals.py --input /Volumes/data/CodexTemp/T2/unfair-01-polynomials/reproduced-raw.json --roots-input H8-a2b3-inputs.json --balls H8-a2b3-root-balls.json --output /Volumes/data/CodexTemp/T2/unfair-01-polynomials/reproduced-certificate.json
python3 run_parent_checker_negatives.py --package . --checker verify_H8_exact.py --output /Volumes/data/CodexTemp/T2/unfair-01-polynomials/reproduced-negative-results.json
```

原始完整有理区间中间件约 15 MiB，不归档。紧凑证书重建 SHA-256 必须等于 `1507a8cc69171d3cf565776c2726b640972107d58cd152f3bf2a7063e408a919`。14 项深层篡改回归实际全部被拒绝；输入/根盘篡改时同步更新外层哈希，以检查语义而非仅检查文件哈希。
