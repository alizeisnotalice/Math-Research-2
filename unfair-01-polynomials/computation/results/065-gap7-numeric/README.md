# 第 7 项缺口的 Metal 与 MPS 数值探路

本机 PyTorch 2.15.0 MPS 完成 `(48,56)` 的 12,288 条并行轨迹与 36,000 步，保存的 288 组 Float32 样本按精确二进制有理数重算。最小保留卷积残差为 14807277/67108864。该次计算进程已退出。

MacBook Apple M5 Pro 的原生 Metal 完成 35 对次数，每对 32,768 条轨迹、32,000 步。保存样本由本机 Python 3.13.18 重算 560 组卷积。原生实现见[第 63 轮程序](../063-gap7-numeric/README.md)。整批监督进程已退出。下一批 40 对不同次数正运行。两类结果均用于选择搜索入口。

```sh
python -B verify_metal_batch.py mac-sweep35 --scratch /Volumes/data/CodexTemp/T2/a30-metal-samples-new --out /Volumes/data/CodexTemp/T2/a30-metal-verification-new.json
python -B ../060-gap7-numeric/verify_saved_samples.py --out /Volumes/data/CodexTemp/T2/a30-mps-verification-new.json .
```
