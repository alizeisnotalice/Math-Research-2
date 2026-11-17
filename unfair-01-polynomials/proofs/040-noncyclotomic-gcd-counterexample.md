# 第 40 轮：二元直接和边的非循环共同因子反例

**状态：严格局部命题否定。** 这不是第 39 轮的已采用结论，也不构成原 0--1 公平猜想的
反例。它否定的仅是下述可能用来处理共同因子分支的辅助断言：

> 若 \(U,V\in\{0,1\}[x]\)、\(U(0)=V(0)=1\)，且 \(UV\) 仍二元，则
> \(\gcd_{\mathbb Q[x]}(U,V)\) 只含循环因子。

令

\[
 \begin{aligned}
 U&=1+x^3+x^7+x^{11}+x^{14},\\
 V&=1+x+x^{13}+x^{23}+x^{33}+x^{45}+x^{46},\\
 \ell&=x^{10}-x^9+x^7-x^6+x^5-x^4+x^3-x+1.
 \end{aligned}
\]

## 自包含的精确证书

令 \(A=\{0,3,7,11,14\}\) 与
\(B=\{0,1,13,23,33,45,46\}\)。\(A+B\) 有恰好 \(5\cdot7=35\) 个不同元素，故
\(UV\) 的每个系数皆在 \(\{0,1\}\)。精确 \(\mathbb Z[x]\) 恒等式为

\[
 U=\ell\Phi_5,
 \qquad
 V=\ell\Phi_6\Phi_7\Phi_{30}\Phi_{66},
 \qquad
 \gcd(U,V)=\ell.
\]

并且 \(\ell(0)=1\)、\(\ell(-1)=-1\)。介值定理给出 \((-1,0)\) 内的实根；循环因子的
根全在单位圆上，故这个 gcd 必含非循环因子。`counterexample_exact.py` 用 SymPy 的精确
整数运算重新核验支撑和的单射性、全部乘积系数、两条分解、gcd 和端点符号。

## 论文线索、实际访问与定位

这个具体的整体 Newman 多项式（\(UV\)，次数 60）来自下列论文的显示式 (7)：

M. J. Mossinghoff, *Polynomials with Restricted Coefficients and Prescribed Noncyclotomic
Factors*, **LMS Journal of Computation and Mathematics** 6 (2003), 314--325,
DOI [10.1112/S1461157000000474](https://doi.org/10.1112/S1461157000000474)。出版商的
[Cambridge 文章页](https://www.cambridge.org/core/journals/lms-journal-of-computation-and-mathematics/article/polynomials-with-restricted-coefficients-and-prescribed-noncyclotomic-factors/3DC031DC40018F9509F93F744CACCC3A)
与[公开 PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/3DC031DC40018F9509F93F744CACCC3A/S1461157000000474a.pdf/polynomials-with-restricted-coefficients-and-prescribed-noncyclotomic-factors.pdf)
均已实际访问。

本轮实际读取 PDF 首页的元数据（作者、题名、期刊为 *LMS Journal of Computation and
Mathematics*）和出版商文章页的 citation 元数据（卷 6、314--325、DOI）；并读到第 2 节
显示式 (7)（印刷页 315；其展开式和 \(\ell^2\Phi_5\Phi_6\Phi_7\Phi_{30}\Phi_{66}\)
分解）。直接 DOI 短链接在浏览器工具中未能解析，但出版商的 DOI 元数据和链接已核。
随后下载同一 PDF 到本目录，仅作临时来源核对。下载文件 SHA-256 是

```text
d578f4a4ad1cd252db412ac568bed20382819278c67ed6f0d92982ac9dbffe0f
```

其 PDF 数字页数为 12，印刷页锚点为 314--325，均已实际核验。论文仅提供发现线索；
上节的支撑、分解和 gcd 都由本地精确脚本独立重算，因而不依赖对来源公式的信任。

## 影响范围

二元边与其直接和乘积本身不足以消除 Lehmer 型非循环共同因子。共同因子路线若要继续，
必须额外使用完整交叉模型的切片兼容性、根选择或正性约束。
