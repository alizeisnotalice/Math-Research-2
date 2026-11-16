# r021 最小自包含归档

本包保留最终调度器、严格汇总器、负测试、独立权重后处理器，以及两次实际小规模
程序核验。它不含 GPU binary、编译缓存、历史 v1--v3 重复文件或大尺寸生产数据。

## 内容

- `source/`：最终 Python 源；汇总器要求并重放相邻 manifest。
- `data/n3/`：实际 (n=3,P=17,R=64,L=4) v4 raw、manifest、严格汇总和负测试。
- `data/n32/`：root 实际 (n=32,P=17,R=8,L=4) 的 raw、manifest 和严格汇总。
- `inputs/tree-audit-n3.json`：原 n=3 整数精确输入，供有理求和核验。

真实 v4 n=3 命令：

```sh
python3 -B source/random_tau_design.py --n 3 --reps 64 --P 17 --moveL 4 --method strat \
  --tau-domain 5210000000000321 --trajectory-domain a7210000000000342 \
  --gpu-binary /absolute/path/to/frozen/n32_extension \
  --out-jsonl data/n3/raw.jsonl --n3-exact inputs/tree-audit-n3.json
python3 -B source/aggregate_random_tau.py data/n3/raw.jsonl --n 3
python3 -B source/test_aggregate_guards.py data/n3/raw.jsonl
```

归档 raw 已由上述实际命令产生；为避免覆盖它，复现 GPU 时应首先使用新输出路径。

`data/n32` 的 8 个非死观测均发生 Double 下溢，但严格 log 汇总仍可运行。其权重
ESS 约为 1 意味着样本几乎被一个重要性权重支配；所报微小 delta-SE 因此不表示有
效精度，更不能作为原条件律概率、渐近结论或反例证据。
