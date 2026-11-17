# 有限 first-merge RB 校准

包内保存 n=3/n=7 的原始有限记录、派生源、最小 n=3 trace、精确 top2-D4 基线与
输入拒绝测试。聚合器强制检查 replicate ID、run_end、seed、每个 RB 记录、群组恒等式、
派生源 SHA，并让 `n3.json`/`n7.json` 与 `summary-n3-n7.json` 的 (A_2,B_{top2}\)
完全一致。错误的 (B\) 基线与缺少 RB record 均被负例拒绝。

n=7 的有限估计仍有 (B\) 相对误差约 .00610488、(A_2\) 相对误差约 .00580131，
且 (Q_{mRB}\approx.3447602283\) 而精确 (Q\approx.3444183148\)。故 `audit_passed`
只表示输入/源码/精确基线的一致，绝非数值精度、置信区间、混合或原题结论。
