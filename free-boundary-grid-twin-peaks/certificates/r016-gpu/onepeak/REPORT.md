# r016归档：一峰 normalizer 大 n 数值尺度批

运行目录内部名为r014，实际于本课题第16轮归档；不是第14轮已交付的n8/n9数据。

本目录是对 \(Z_1(n)=A_1(n)/(n^2)!\) 的有限 SMC 计算。理论 FK 恒等式、
一峰与连通前缀的等价、以及理想随机数下的一阶无偏范围见 `THEORY.md`。
它不输出一峰的 \(U\)：两峰的第二根/距离矩在这里没有定义意义。它也不是
关于 \(R_n=A_2/A_1\)、远峰概率或 \(n\to\infty\) 的严格结论。

## 冻结核和锚定

去除Metal字符串首尾换行后的SHA-256是
`632e9990f8f2721efcd02f7f8b4ec1c755ce5626589697c7caccbd10fdd8167e`，与
r013的同样截取内容逐字节相同。保留首尾换行的标准提取SHA-256为
`7474e22aee24049eb8692db789f32dff7b46e5c1c9d31f569622a7278c47da39`；主Agent独立核对标准提取内容。host 固定 `tau=N+1`，故 rank 1 在全图、其后
各秩在 boundary；full single-command-buffer integer-stratified mode 的 parent
目标使用 domain 1，选 parent 后的 mutation 顶点使用独立 domain 2，且每重复
在 autoreleasepool 内执行。

按已归档精确整数

\[
A_1(7)=134267527713152677457983522684235911142268448963584,
\qquad A_1(7)/49!=2.2073241970861836\times10^{-13}.
\]

r014 的 n=7 小链路校准（P=8192、R=16、domain `8300…`）均值是
\(2.2250933725939906\times10^{-13}\)，相对差 \(+0.00805010\)。它只核验
新目录、参数和归一化链路；主生产使用下表指定的 P/R。

## 实际大粒子批

每行均为 P=262144、R=256、单一 \(\tau=N+1\)，无死端。均值、样本方差与
描述性 SE 均由该文件的256个保存 \(Z\) 重复直接计算。两个新 domain 的差异
只是有限 PRNG 数值尺度检查，并不证明重复独立或建立严格置信区间。

| n | domain | \(\bar Z_1\) | 样本方差 | 描述性 SE | wall / GPU时间和 (s) |
|---:|---|---:|---:|---:|---:|
| 10 | `8300…` | 9.4264807888e-31 | 3.8279567135e-63 | 3.8669052112e-33 | 176.442 / 175.018 |
| 10 | `8400…` | 9.4408663772e-31 | 3.4288866420e-63 | 3.6597934976e-33 | 161.411 / 159.887 |
| 12 | `8300…` | 1.6341870642e-46 | 4.3440097010e-93 | 4.1193188629e-48 | 274.382 / 272.387 |
| 12 | `8400…` | 1.6068093804e-46 | 2.3019368576e-93 | 2.9986565075e-48 | 239.796 / 238.009 |

GPU 时间是每个 command buffer hardware timestamp 的跨重复和；同步 submit+wait
时间和、GPU 时间和和全进程 wall 属于不同聚合层次，不能相加。每一 raw 保存所有
seed、每重复 Z、死端和计时；四份均为 258 行（metadata、256 repeats、run_end）
且以退出码0完成。进程回收已核验。

## 复现与限制

```sh
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules
swiftc -O -framework Metal source/onepeak_smc.swift -o onepeak_smc
./onepeak_smc 12 262144 256 8300000000000000 raw/n12-p262144-r256-d8300.jsonl
python3 source/summarize_onepeak.py 'raw/n1[02]-p262144-r256-d*.jsonl'
```

归档文件hash见课题根目录 `manifest-r016.json`。有限 PRNG、有限粒子和有限重复意味着上表仅是
数值证据；未保存完整标号，不能用这些样本替代独立几何验证。与另一 worker 的双峰
normalizer 比较 (R_n) 时，还必须保留两个估计的 proposal、seed、方差与偏差限制，
不能直接把表中小数视作严格比值。
