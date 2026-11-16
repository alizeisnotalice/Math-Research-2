# 给定16339个因子的有理不可约性证书完成

本包补齐第36轮剩余4项。原16339项中的13149项见第33轮，3186项见第36轮，本包4项完成总目录。
[完整论证](proof.md)给出九次核的有限域排除、平方根扩张、互反迹多项式与非平方判别式；有限域约化用迹矩阵分母解释，不依赖整数环软件。
verify_four_qq.py只用整数运算，父版本补查源恰为4项及两个平方复合恒等式。父实际复跑0.0031秒，五类语义负测通过；Agent原回执独立保留，不混同其源哈希。

```sh
/Users/zzh/.venvs/main/bin/python verify_four_qq.py --bound ../036-factor-degree-sieve/bound-verified.json --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/four-qq-replay.json
```

仅给定有理因子目录认证完成，不是全部实因子布局认证或原猜想证明。
