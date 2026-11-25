# 欧拉常数的无理性
研究对象为实数 γ = lim(H_n − log n)，其中 log 为自然对数。目标是对所有互素整数 a,b（b>0）排除 γ=a/b。首轮复现 Sondow 判据与 Aptekarev 逼近，追踪清分母后的严格余项。

## 研究入口

- [第 1 轮中文报告](rounds/001.md)给出两条证明路线、有限结果、证明主张核查与后续引理。
- [27 项阅读索引](references/index.md)、[BibTeX](references/references.bib)、[查询记录](references/search-log.md)。
- [严格区间算法](proofs/intervals.md)、[消去对数的整数关系](proofs/blocks.md)、[当前状态](STATE.md)。

## 重放

在本目录下执行。精确证书仅需 Python 标准库。本轮在 macOS 上用 CPython 3.13.18、mpmath 1.3.0、SymPy 1.14.0 运行；精确证书另以 CPython 3.15.0 重放。全程确定性。数值积分在 22、32 位十进制精度分别运行，允许绝对差 10^(7−精度)；最小二乘用 80 位。仅数值脚本依赖 mpmath，块实验依赖 SymPy。

~~~sh
python3 -B computation/replay.py --check computation/certificate.json
python3 -B computation/intervals.py 2 -1
python3 -B computation/intervals.py 0 0 --log 2 4 --mode floor
python3 -B computation/approximants.py 12
python3 -B computation/sondow.py 12
/Users/zzh/.venvs/main/bin/python3 -B computation/numerical.py
/Users/zzh/.venvs/main/bin/python3 -B computation/blocks.py
~~~

intervals.py 的两个位置参数分别是 γ 的系数和整数常数；可重复使用 --log K COEFF。返回有理端点、精度和严格符号或 floor。approximants.py 输出 n=0 至给定上限的 JSON，error_interval 是带符号的 γ−p/q，form 是 qγ−p。sondow.py 输出单个指标的素因子指数、log S 的区间与 floor、dI 的区间、排除的分母所属的因子集。scaled_threshold_difference 是 \(2^n(\{\log S_n\}-2^{-n})\)。推荐先重放 n=1 和 n=12。

[certificate.json](computation/certificate.json) 保存有限结论的精确端点。[numerical.txt](computation/numerical.txt) 保存积分与拟合输出。[blocks.json](computation/blocks.json) 保存块实验的精确坐标。需要本机重定向额外输出时使用 /Volumes/data/CodexTemp/T2/euler-gamma-irrationality/。
