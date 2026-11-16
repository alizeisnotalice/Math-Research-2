# 首分数位置的完整低端前缀与MPS搜索原型

PROOF.md证明低端长度d的精确充要条件及两个全局必要筛，均不保证全局延拓。d1..9完整账本已由独立四状态枚举 verify_prefix_ledger.py 重算；主生成器用三状态并非验证器的覆盖来源。

```sh
python3 verify_prefix_ledger.py --input fractional-prefixes-d1-d9.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/prefix-recheck.json
python3 run_prefix_batch.py --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/prefix-new-run --m 32 --n 72 --d 7 --starts 16384 --steps 12000 --restart 2000 --seed 2026111605 --margin 0.001 --seconds 14400
```

GPU脚本需本机PyTorch2.15.0、NumPy与MPS。完整d7批次属于下一轮，不计本包的完成覆盖。driver-smoke为实际m12/n18/d5、512初值、200步、100步重启测试；重启前旧16样本与最终16样本已保存，前缀从系数独立重建。驱动退出0，父子均reaped/不存在。

每500步及结束保存CPU state.pt，包括优化器与张量；可再生state.pt不入库，未保存完整随机数状态，不声称逐比特恢复。重启前保存旧组，浮点互补和仍须精确检查。检查点和低残差不构成反例或有限排除。prefix-gpu-load-samples.json是后继真实d7批次连续4次设备100%采样，仅代表采样窗口。
