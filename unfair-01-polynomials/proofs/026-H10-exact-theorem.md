# 固定目标 \(H_{10}(x^2)H_{10}(x^3)\) 的严格非负分解定理

令 \(H_{10}(x)=1+x+x^2+x^9+x^{10}\)，\(C=H_{10}(x^2)H_{10}(x^3)\)。

## 定理

对任意 \(A,B\in\mathbb R_{\ge0}[x]\)，若 \(A(0)=B(0)=1\) 且 \(AB=C\)，
则
\[
 \{A,B\}=\{1,C\}\quad\hbox{或}\quad\{H_{10}(x^2),H_{10}(x^3)\}.
\]
这是固定的次数 50 二元目标的结论，**不**解决一般 0--1 多项式猜想。

## 证明链

模 \(7,73,19\) 的精确平方自由因子型 \([10]\)、\([1,9]\)、\([2,3,5]\)
认证 \(\operatorname{Gal}(H_{10})=S_{10}\)：传递性、点稳定子中的 9-循环和
换位依次给出二重传递与全对称群。\(H_{10}\) 非自反，故第 18 轮根幂积平衡
引理适用。

最终证书以有理 Rouché 不等式和分离圆盘认证五个正虚部母根及其共轭、以及每一
母根的两个 \(x^2\)-提升、三个 \(x^3\)-提升及其共轭：共 50 个互异目标根。
若 \(AB=C\) 如上，非负性和常数项一使每个因子系数不超过 C 的对应二元系数，
最高系数乘积为一，故 \(A,B\) 首一。于是 \(A\) 选取 25 个实二次根群的一个
子集，\(B\) 取补集。

把十个母根分为五个共轭块。第 r 块选取 \(u_r\in\{0,1,2\}\) 个 a-lift 群及
\(v_r\in\{0,1,2,3\}\) 个 b-lift 群。将 A 的根积取六次幂并用根幂积平衡，
得到共同 E 满足 \(3u_r+2v_r=E\)（r=1,...,5）。这只是实际非负因子必须满足的
必要条件。按补集取代表，非恒定布局数为
\[
 3^5+2^5+3^5+6^5+2^4=243+32+243+7776+16=8310.
\]
完整量词及终端补集约化见 [覆盖说明](026-H10-root-layout-coverage.md)。

对每个根圆盘，checker 用有理闭盒包住
\((t-z)(t-\bar z)=t^2-2\Re(z)t+|z|^2\) 的系数；它以 Fraction 端点乘法重算
每个布局或补因子的指定系数。8310 个独立生成布局中，8309 项的认证上端严格
小于零，故相应首一根积有负系数，不可能对应非负因子。最小认证余量为
\[
 257453695/1099511627776.
\]
唯一未排除的 \((E,\mathrm{code})=(6,0)\) 精确选取全部 20 个 a-lift 根，
补集选取全部 30 个 b-lift 根；经根数与首一性识别为
\(H_{10}(x^2)\)、\(H_{10}(x^3)\)。空集/全集给常数分解，定理得证。

## 固定证书与复现

独立 checker：`verify_H10_exact.py`，SHA-256
`57354c1804d9dbb7ce631afa92047a2d01bbbf19468432c3bb9ba03a1113f7a1`。
已通过的输入、balls、紧凑证书、母式模证哈希依次为
`15bad8ecbb194b3f12b2b550adcd0d1a04b29f457a6b21563eebfbacae264332`、
`2019946b5d7637e4ec799b66ef76b1f3462eec5502e33876f5b221f8e750e8d3`、
`b3c4fae2a9b72fd8951441f9065bf05094ac23694c149a1151ef50c949cdb3ae`、
`a6559a16f8ca29459f452288b1b62e004a9b809ceb9f2bb543b0c6e0e591d97f`。

从仓库根运行，生成器和严格检查器的同名 helper 分别放在 generation/ 与包根：

```sh
cd unfair-01-polynomials/computation/results/026-H10-exact
python3 verify_H10_exact.py --input H10-a2b3-inputs.json --balls H10-a2b3-root-balls.json --certificate H10-a2b3-certificate.json --galois H10-mother-certificate.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H10-recheck.json
```

白名单：本稿、覆盖说明、独立 checker 与两个精简 helper、
`H10-a2b3-independent-verification.json`，以及
`test_H10_checker_guards.py` 和其 v4 收据。后者实际拒绝重复 negative/unknown、
零半径、family/slot 错配、非正虚部/重复母根、记录遗漏、错 digits、伪正上端、
伪标准项和错误有限域因子；它们不是仅靠 hash 不匹配拒绝。
