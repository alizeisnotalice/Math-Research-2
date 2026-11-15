# r007 high-beta GPU handoff

This is a construction search, not a conditional-probability experiment.

## What survives

`r007-metal.jsonl` was overwritten by the v2 run and retains 16 complete
samples only: eight at n=24 and eight at n=32.  The earlier n=8 smoke and n=12
v1 batches have only their printed summaries in the chat/log history; their
full labels are not recoverable and are not claimed here.

The independent audit `audit_r007.py` checks every retained label is exactly
a permutation of 1..n², recomputes two peaks, distance, integer M,b,t, beta,
first merger and steepest-ascent basin interface.  Its actual output is
`audit-v2.json` (16/16 passed).

| n | distance min–max (mean) | beta min–max (mean) | interface edges min–max (mean) |
|---|---|---|---|
| 24 | 12–29 (19.375) | .949477–.982578 (.968859) | 35–47 (40.375) |
| 32 | 25–62 (38.375) | .941292–.977495 (.963546) | 16–58 (43.875) |

Thus accepted swaps do not force the seed to remain literally one thin path:
some retained samples have substantially larger gradient-basin interfaces.
They nevertheless remain a finite high-beta family, with no count relative to
A2 and hence no normalized lower bound.

## Versions and timing

* `r007.swift`, SHA-256 `38844497...131d9113e`: v1; only submit-to-complete
  wall time was accumulated.  Its historical field named `gpu_seconds` must
  be read as wall time, never hardware execution time.
* `r007_v2.swift`, SHA-256 `a6bd1cdf...c6f5b8590`: adds
  `cb.gpuEndTime-cb.gpuStartTime`, after SDK verification that both properties
  exist.  `r007-v2-timing.txt`, SHA-256 `08cf5549...04ab8b21`, records each
  dispatch separately.
* For v2 n=24+n=32, 32 dispatches have wall sum 29.5487997532 s and Metal
  hardware-execution sum 29.5421959999 s; every hardware duration was >0.
  Both n values used 2,048 chains x 16 rounds x 2,048 proposals = 67,108,864
  proposals each.  n=24 accepted 15,603,564 swaps and best beta .9877408056;
  n=32 accepted 12,616,491 and best beta .9872673849.

No binary is included.  All GPU processes were checked exited after audit.

## 主Agent逐文件核验与更正

此前handoff目录只有README和audit-v2.json；主Agent已从原运行目录实际复制本报告涉及的两版Swift源、现存16个标签、v2逐dispatch计时、CPU审计及种子核验脚本，并生成SHA256SUMS。

v1/v2源码真正使用的随机种子基数为0xA0227000，即2686611456，逐链加链号。历史summary中的硬编码seed=2686607360错误，差4096；保留历史JSON，复现以实际源码种子为准。

Swift每轮在CPU重算全部链的K和距离并计算β；并未把GPU输出buffer中的β与CPU β逐项断言相等。因此这里不声称完成了GPU/CPU β数值一致性检查。保存样本的β及M,b,t均由独立CPU重新计算，样本有效性不依赖GPU评分数值。

v2记录gpuEndTime-gpuStartTime，但宿主未显式检查command buffer最终status；这些运行是探索证据，不能作为认证采样。后续采样程序应检查status及所有跨设备指标。源码对现有输出文件有覆盖行为；复现须先创建独立输出目录和空的r007-metal.jsonl，不能在归档目录运行。

主Agent还补入了保存distance及整数M,b,t与独立重算值逐项相等的断言，并实际重跑16个样本；与原audit-v2.json全部一致。
