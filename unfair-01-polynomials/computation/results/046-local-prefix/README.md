# 第 46 轮：N=5 局部前缀反例包

这是对辅助入口“由 \(A_5=B_5=0\) 和前五个二元卷积行排除 \((r,t)=(1,5)\)”的精确反例。
它不是原均匀和问题反例。

精确数据为

\[
A=1+\frac9{10}x+\frac{181}{1000}x^3+\frac{2421}{40000}x^4,
\]

\[
B=1+\frac1{10}x+\frac{91}{100}x^2+\frac{7371}{8000}x^4.
\]

两者常数项一、逐系数非负，且 \(A_1,B_1>0\)、\(A_5=B_5=0\)。精确卷积为

\[
[x^k]AB=1\ (0\le k\le5),\qquad
[x^6]AB=\frac{220311}{4000000}\in(0,1).
\]

因此它仅说明这个有限局部推理入口不足。第六行已不是二元，且严格分数有理系数不符合
原问题中首一整数二元目标因子的代数整数系数条件。

`verify_r1t5_local_prefix.py` 与 `verify_r1t5_local_prefix.wl` 是已实际运行的两份精确复算；
对应 stdout 回执一并保留。复现命令：

```
python3 verify_r1t5_local_prefix.py
wolframscript -file verify_r1t5_local_prefix.wl
```
