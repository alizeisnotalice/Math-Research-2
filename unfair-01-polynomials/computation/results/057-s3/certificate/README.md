# 第 57 轮八次位置 3 的完整证书

完整链覆盖 13,209 个实二次核心和 1,915,327 个整尾组。1,523,605,945 个必要母式以整数区间表示。双嵌入正轴、Cayley 线性界与严格根环留下 869,396 个母式。正向树排除 869,395 个，剩余一个在反向第 26 层排除。未知为零。

- `data/` 保留排序的所有根筛通过项及双树逐母式账本，压缩分片逐字节可重建。每条树记录含母式、排除方向、层数、首次空层、最大活状态数及逐层记录的哈希。
- `original/` 保留实际根筛计数、边界回归、全量双树与第三套预算树回执。第三预算树的六个整数上限触发使用任意精度回退，并逐母式核对排除方向、深度和最大活数。
- 正轴和 Cayley 按尾区间相交。无需保存 15 亿个母式逐点文件。鲜跑重造完整尾组、根通过项、双树与第三预算树，并比较输入及确定性输出的全部承诺。

使用 Python ≥3.11、Clang C++20 和 GMP。输出目录须为证书外的新路径。纯归档模式验证字节与覆盖。实际全重跑用第二条命令。

```sh
python3 -B modules/verify.py
python3 -B modules/verify.py --workers 8 \
  --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay057-new \
  --gmp-prefix /opt/homebrew \
  --clang-module-cache /Volumes/data/CodexTemp/cache/shared-clang-modules
```

根筛编译禁用 fast-math 和 FMA contraction。程序检查 IEEE binary64 和最近舍入模式，并用严格边界及次正规数回归验证当前环境。所有区间含糊项走 GMP。第 57 轮独立定点整数区间根筛见主线对应证据。
