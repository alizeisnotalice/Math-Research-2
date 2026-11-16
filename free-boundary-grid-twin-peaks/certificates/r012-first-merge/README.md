# 首次合并空间矩和根深度的整数账本

本轮新增 n=3,4,5 的完整首次合并层表；此前的 A2、Q 精确值是独立基线，不是本轮新发现。

对于非空集合 C，a(C) 是内部唯一峰的排序数；Mx、My、M2 是该根的坐标及平方范数总和。以单点 (1,x,y,x²+y²) 为初值，删除最后出现顶点 v 时要求 v 与 C\{v} 有邻接，递推四个量。h(V)=1，h(S)=Σ_{v∈∂S}h(S∪{v})。

枚举连通 S 中使 S\{w}=A⊔B 恰为两个分量的 w，A 固定为含最小编号的分量，不另乘2。计重为

\[
W=h(S)\binom{|S|-1}{|A|}a(A)a(B),
\]
\[
D=h(S)\binom{|S|-1}{|A|}
\{a(B)M_2(A)+a(A)M_2(B)-2[M_x(A)M_x(B)+M_y(A)M_y(B)]\}.
\]

每一双峰降序标号有唯一首次合并状态；两分量内部都是单峰，合并前的标签交错任意，合并后不能再出生。因此 ΣW=A2、ΣD=n²A2 Q。这是此前唯一分解的完整有限执行。

## 算术范围

n≤5，s=|C| 时 a≤s!，坐标≤4、平方范数≤32，故 Mx,My≤4s!、M2≤32s!，h(S)≤(N−|S|)!。每个状态的 W≤(N−s)!(s−1)!；D 的正项≤64倍此量。整个合法标号总平方距离≤32·25!；根深度平方双根和≤2·64·25!。所有这些以及递推的非负中间和均远小于2^128。choose64 的 n≤24，其逐步乘法也在64位内。代码没有用128位平方 Mx。所有减法先核对非负。

n5共2,301,812个连通S、13,665,496个合并事件、396个非零(合并秩,小分支大小,跨支最大曼哈顿跨度)层。跨度是 max_{u∈A,v∈B}d1(u,v)，不等同于实际根距。20个最大W状态附原始矩；等权状态取舍可随线程分块变化，不把其列表次序作为认证对象。

## 根深度证书

令 δ_A(p)=d_G(p,V\A)，使用原完整方格距离。两支不邻接，所以沿任意 p 到 q 的路径，离开A与进入B之间至少有一个位于两者之外的点。于是 ρ≥δ_A(p)+δ_B(q)。H(d) 为两根深度的合并计数，严格有

\[
\sum_d H(d)=2A_2,\quad
Q_n\ge\frac{\sum_d d^2H(d)}{2n^2A_2},\quad
\mathbb P_2(\rho\ge r+1)\ge\frac{\sum_{d\ge r}H(d)}{2A_2}.
\]

最后一式因每标号最多贡献两个深根、而另一根深度至少1。指定根计数a_r只改单点初值；逐D4根轨道计算并乘轨道大小。n5为6轨道，全部深度计数和精确等于2A2。

| n | 深度给出的Q下界 | 下界占实际Q |
|---:|---:|---:|
|3|2699/17976 ≈ .15014464|.30670|
|4|4569927227/48706198380 ≈ .09382640|.22190|
|5|70701713638102136309/1103634606462962231250 ≈ .06406261|.16726|

这些是有限严格结果。下降趋势不能判定渐近方向。

## 实际验证和资源

主Agent重新编译空间矩程序并运行n5、12线程，ΣW和ΣD及396层整数逐项等于导入包；墙钟2.51s、CPU user10.83s、sys0.57s、最大RSS2,755,870,720 bytes、零swap，进程结束。源运行n5墙钟2.04s、RSS2,755,837,952 bytes；两次不要合并为一次测量。

根深度源运行n5用12线程、3块128位数组，墙钟5.67s、user23.75s、RSS1,680,769,024 bytes、零swap，进程已退出。n3独立9!枚举直接提取首次合并和两根深度，与DP逐项相同；主Agent再运行该审计及全部整数账本核验。软件为本机Apple clang21/C++20及Python3；环境细项见本轮总清单。

从本目录复现，编译产物和输出可改到CodexTemp：

```sh
c++ -O3 -std=c++20 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang first_merge_spatial_dp.cpp -o /Volumes/data/CodexTemp/T2/first_merge_spatial_dp
/Volumes/data/CodexTemp/T2/first_merge_spatial_dp 5 first-merge-n5.json 12
python3 -B verify_first_merge_spatial.py
cd depth
c++ -O3 -std=c++20 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang first_merge_root_depth.cpp -o /Volumes/data/CodexTemp/T2/first_merge_root_depth
/Volumes/data/CodexTemp/T2/first_merge_root_depth 5 root-depth-n5.json 12
python3 -B verify_root_depth.py
python3 -B audit_root_depth_n3.py
```

`source-index.json` 保存本轮导入原文件路径与SHA256；不归档二进制、全子集数组或缓存。
