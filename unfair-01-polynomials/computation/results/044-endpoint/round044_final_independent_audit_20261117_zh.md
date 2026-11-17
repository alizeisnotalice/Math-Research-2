# 第 44 轮最终独立审查（只读父材料）

审查对象为同步证明 `044-fractional-endpoints-and-quadratic-core.md`、
`044-local-prefix-counterexample.md`，以及父临时目录的 `verify_local_N8.py` 与
`certificate.json`。本审查未改写任何父文件。

## 结论

通过，未发现会扩大结论范围或破坏论证的缺口。

1. 双端分数指标只使用整数前缀、非负性和同时反转；没有单独反转一个因子。
2. \((r,t)=(1,4)\) 的上界正确。式 (7) 中的三次因子确实正：在 \([0,1]\)，
   \(5u^3-3u^2\ge-4/25\)，故再加 \(3u+3\) 严格正。
3. 当前主稿第 118 行已显式写明 G、Q、P 首一且常数项一，足以支撑高端三角展开；归约段
   未假设 \(G/Q\) 多项式或正，也未假设 P 或共轭 P 系数非负。
4. N=8 例只满足截断前缀：它有 \(C_0=\cdots=C_8=1\)，但 \(0<C_9<1\)。因此它只否定
   全 N 的局部前缀引理，绝不是原均匀和问题反例；父稿范围标注准确。

## 实际独立复算

本目录 `independent_verify_local_N8.py` 只读父 `certificate.json`，以标准库
`fractions.Fraction` 重算最小多项式、隔离区间、全部 A 系数区间、商环内
\(C_0=\cdots=C_8=1\)，以及 \(0<C_9=tA_5<1\)。实际输出：

```
PASS: exact Fraction audit: root isolation, positivity, C[0..8]=1, C9 in (0,1)
```

当前主稿数学显示环境静态检查为 `\[` 25 个、`\]` 25 个，未见缺失闭合符号。

## 当前对象 SHA-256

* core proof: `cc34cb668103824b9768a8004dd884f660080fd6b5293115c8fa9e614cd6129b`
* local-prefix proof: `10583ddb87e393e2a5b8e7717cc948d636a7328ca2f8fbff00f9030a21bb4011`
* parent Fraction verifier: `9210e7650be894ce95973621341926de2749dd958e579d32beb087151e2a634c`
* parent certificate: `187545af24ae88d57e7518f81ffbb45d61ed8c01d4367b5233cfb226774222a2`
* independent verifier: `aa06ffbfbfd32c44bb6247e04de727a38114b8ecb8d1e141ee82d59d5b839514`
