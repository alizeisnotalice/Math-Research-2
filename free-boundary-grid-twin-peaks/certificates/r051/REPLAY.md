# 第51轮复现

实际环境为 ARM64 macOS、Python 3.15.0、Apple clang 21.0.0 和现有 Swift Metal 程序。在 data 卷选择空目录作为 SCRATCH，编译共享缓存位于 /Volumes/data/CodexTemp/cache。下列命令从本证书目录运行。

## 低支删除、换根与依赖间隔

表的列为 `i b rho ell h_star count`，\(\ell\) 是低峰到 \(B\cap N(w)\) 的支内最短距离。固定 \(B\) 时，令 \(S\) 表示已见且在 \(G\setminus B\) 连通的高支前缀，初态 \(S=A\cup\{w\}\)。令 \(C_B(S)\) 为 \(B\) 的外邻中尚未见、且在 \(G\setminus B\) 没有已见邻点者，\(h_G(U)\) 为全图从 \(U\) 延伸的连通后缀数。首个依赖间隔计数满足
\[
 f^B_1(S)=\sum_{v\in C_B(S)}h_G(B\cup S\cup\{v\}),\qquad
 f^B_{j+1}(S)=
 \sum_{v\in N_{G\setminus B}(S)\setminus S}f^B_j(S\cup\{v\}).
\]
若 \(C_B(S)\) 为空，随后保持为空。有限间隔数之和从 \(h_G(B\cup S)\) 扣除，即得到完全可删除的计数，输出哨兵 \(h_*=N+1-i\)。程序用两层滚动数组并对 \(B\) 取正方形对称轨道。

```sh
clang++ -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang \
  -std=c++20 -O3 -Wall -Wextra -Werror essential/essential.cpp -o "$SCRATCH/essential"
"$SCRATCH/essential" 3 4 1 8 | cmp - essential/n3-wide.tsv
"$SCRATCH/essential" 4 4 1 8 | cmp - essential/n4-wide.tsv
"$SCRATCH/essential" 5 4 1 8 | cmp - essential/n5-wide.tsv
python3 -B essential/verify.py | cmp - essential/verification.json
python3 -B strip/audit_wide_strip.py | cmp - strip/audit_wide_strip.stdout
python3 -B reroot/audit.py | cmp - reroot/audit.stdout
```

后三个位置参数是线程数、低支大小下限、上限。五阶生产分四段，大小1和3–5在本机，大小2和6–8在 MacBook。1、2使用原冻结源码 essential_b2_frozen.cpp，3–8使用最终源码。远端包装用 run_b2.py、run_wide.py。两份完成回执中的原始 TSV 哈希可由合并表逐区间筛出。验证器对照三阶所有 \(9!\) 个词，第50轮删除计数及第49轮首接触分布。短路径逆像独立检查三阶217个连通真子集及每种编码的80,880个接触与词对。

## 十二阶严格有限界

[三矩形与连接秩证书](n12/N12_T_BOUND.md)。将 rect_protected.cpp、count_rectangles.py、numerator_n12.py、pratt-61bit.json 复制到 SCRATCH，设置 A022_RESEARCH_ROOT 为仓库内课题根目录，用原 clang 参数重编译 rect_protected。

```sh
python3 -B count_rectangles.py output-directory --workers 3
python3 -B numerator_n12.py output-directory/summary.json numerator.json
```

已重跑32个模数计数及穷举门槛。比较两份分母的 computed_results、regional_histograms、regional_plans、模数和校准字段。测时、RSS、驱动源码哈希可能不同。已用归档分母 n12/denominator-final.json 再生成分子，与 n12/n12-conditional-final.json 逐字节一致。复现驱动改用可移植路径，原回执保存初次驱动哈希。

## 十二阶本机 Metal

```sh
python3 -B local-metal/verify.py | cmp - local-metal/summary.json
```

第26轮 Swift 源的哈希与完成回执一致。重跑时将源码与 local-metal/run.py 复制至同一个空目录，用 `swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift` 编译为 tau2_merge_rb_calibration_r026，再运行包装脚本。实际参数为 \(n=12\)、131,072 粒子、16 重复、\(\tau=2\)、域 9e51128000000001、分层重采样、128 次移动。短链对照来自第26轮原始日志，读取该日志核对哈希后归档对照汇总。两批种子和重复数不同，作为有限数值敏感性检查。
