# 支配补集晚合并条件：首次合并权重的精确覆盖

日期：2026-11-16（Asia/Shanghai）。本包只做有限精确诊断，不修改主仓或此前冻结的 DP 包。

## 被统计的事件

首次合并的两个规范分支取为 (A,B)（删去桥点后，含最小顶点的分量记为 (A)），并令

\[
C=V\setminus(A\cup B),
\]

所以 (C) 包含首次桥点 (w)。分别测试

\[
P_A:\ \forall c\in C,\ N(c)\cap A\ne\varnothing,
\qquad
P_B:\ \forall c\in C,\ N(c)\cap B\ne\varnothing.
\]

覆盖事件是 (P_A\cup P_B)，并精确分为 only-A、only-B、both。对每个状态仍用首次合并的原始权重 (W) 和根距平方总量 (D)，所以覆盖条件下的
\(Q\) 是 (D/(n^2W)\)，不是抽样频率。

对一个满足侧 (X) 和另一侧 (Y)，令
\(c_0=\#\{c\in C:N(c)\cap Y\ne\varnothing\}\)。每个这样的分割恰给出 `c0` 个桥状态，且逐项实际验证
\(h(S)=(|C|-1)!\)。only-A/only-B 选择相应满足侧；both 只选择规范 (A)，避免一个分割被作为两侧重复。每个 `(special_a, other_b, g, c0)` 层都验证 `state_count/c0` 是整数。

## 覆盖与条件空间矩

| n | 覆盖 W/A₂ | 覆盖内 Q | 覆盖 D/全 D |
|---:|---:|---:|---:|
| 3 | 263/749 ≈ 0.35113485 | 140/263 ≈ 0.53231939 | 0.38181818 |
| 4 | 39851458/193278565 ≈ 0.20618664 | 6180009613/13390089888 ≈ 0.46153608 | 0.22505541 |
| 5 | 3628522440251982214/36787820215432074375 ≈ 0.09863380 | 7914509430702761293/18142612201259911070 ≈ 0.43623869 | 0.11234313 |

n=5 的 only-A、only-B、both 覆盖质量分别为全 (A_2) 的
0.06621809、0.02704008、0.00537562；它们的条件 (Q) 分别约为
0.43208、0.44075、0.46481。

因此该支配补集/晚合并条件在小尺寸中确实承载非零的原始首次合并质量，但覆盖比例从 n=3 的 35.1% 降到 n=5 的 9.86%。这只是有限趋势，既不证明质量趋零，也不反驳存在其他可归一化的晚合并族。

## 引理比值分层（不相加）

每一层输出第13轮旧引理的单分割比例

\[
R=\frac{c_0}{g}\frac{\binom{a+b}{a}}{\binom{N-1}{b-1}}.
\]

n=5 的 247 个分割参数层中，152 层有 (R<1)，95 层有 (R\ge1)。最小实际参数层为 `(a,b,g,c0)=(3,15,7,5)`，其
\(R=10/33649\)。`R<1` 的选侧状态本身覆盖 (A_2) 的约 3.1073%，其余选侧状态约 6.7561%。这些只是按状态的诊断：**没有把各分割的 (R) 相加，也没有由这些有限层作渐近外推。**

无重数分割族账本中，n=5 有 2,598,856 个桥状态、937,133 个规范分割；全部 `h` 阶乘检查和 `state_count/c0` 整除检查通过。

## 独立核验、资源与复现

- n=3 直接枚举全部 (9!\) 标号，筛选 143808 个双峰标号，直接恢复其首次合并 (A,B,C)、覆盖类别和峰距平方；only-A、only-B、both 的 `W,D` 均逐项等于 C++ 账本。
- n=4/n=5 重获既有总 (A_2,D) 的子量界，且类别和、比值分类和、阶乘与整除条件由 `verify_dominating_complement.py` 实际验证。
- n=5 前台 12 线程：wall 2.46 s、user 11.42 s、sys 0.93 s、最大 RSS 2,760,572,928 bytes、0 swaps。运行结束后无残留进程。

复现：

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
c++ -O3 -std=c++20 -pthread \
  -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang \
  dominating_complement.cpp -o dominating_complement
./dominating_complement 3 dominating-n3.json 12
./dominating_complement 4 dominating-n4.json 12
./dominating_complement 5 dominating-n5.json 12
PYTHONDONTWRITEBYTECODE=1 python3 audit_dominating_n3.py
PYTHONDONTWRITEBYTECODE=1 python3 verify_dominating_complement.py
```
