# R059：\(\Psi\) 的 finite-child 与 no-child 机制

## 精确分解

沿既有一峰词记号，令
\[
 \Psi_{\rm fin}:=\sum_{\substack{p\ {\rm eligible}\\b_p\le N}}
       \frac{l_p-1}{b_p-l_p-1},
\qquad
 \Psi_\infty:=\sum_{\substack{p\ {\rm eligible}\\b_p=N+1}}
       \frac{l_p-1}{N-l_p}.
\tag{1}
\]
于是逐词严格有
\[
 \Psi=\Psi_{\rm fin}+\Psi_\infty.
\tag{2}
\]
\(\Psi_{\rm fin}\) 记录至少一个实际 sole-parent child 的父点；R055 的 clean-chain
引理只能覆盖它的一个受限子集。 \(\Psi_\infty\) 记录没有任何 sole-parent child
的父点，不能由 clean-chain 打包替代。

## 两个基本词支付不同的项

逐行词 \(\operatorname{rank}(x,y)=yn+x+1\) 中，所有 \(x,y\ge1\) 的点没有
sole-parent child，故属于 \(\Psi_\infty\)。对
\(y\ge\lfloor n/2\rfloor+2\)、\(1\le x\le n-1\)，令 \(r=n-y\)。有
\[
 l_{(x,y)}-1=(y-1)n+x\ge N/2,
\qquad
 N-l_{(x,y)}\le (r+1)n.
\]
故严格地
\[
 \Psi_\infty\ge
 \frac N2\sum_{r=1}^{\lceil n/2\rceil-2}
 \frac{n-1}{(r+1)n}
 =\frac{N(n-1)}{2n}
 \bigl(H_{\lceil n/2\rceil-1}-1\bigr).
\tag{3}
\]
特别地（\(n\ge8\)），逐行词的 \(\Psi_\infty=\Omega(N\log n)\)。这解释了 R058 的
clean-chain 覆盖失败并不导致 \(\Psi\) 小。

相反，connected-dominating comb 词在其骨架形成后有大量 \(b_p=N+1\) 的补集点，
但这些点的 \(l_p\) 早。事实上有如下通用控制。若某一峰词在秩
\(\tau<N\) 已被前缀支配，则所有 eligible 点的 \(l_p\le\tau\)，从而
\[
 \Psi_\infty
 \le m_\infty\frac{\tau-1}{N-\tau}
 \le N\frac{\tau-1}{N-\tau}.
\tag{4}
\]
这里 \(m_\infty\le N\) 是 no-child eligible 点数；函数
\((x-1)/(N-x)\) 随 \(x\) 增长，故逐项估计即证。

对奇数 \(n\ge5\) 的 comb，前两段大小
\[
 T=\frac{N+2n-1}{2}
\]
已经支配全图，所以 \(\tau\le T\)。代入 (4) 并作保守放松，得到
\[
 \Psi_\infty\le\frac{NT}{N-T}
 =N\frac{n^2+2n-1}{(n-1)^2}
\le\frac{17}{8}N<3N.
\tag{5}
\]
最后一个不等式在 \(n=5\) 取等号，之后单调变强。此前的严格 comb 链调度证明提供
\(\Psi_{\rm fin}=\Omega(n^3)\)。因此两个基本词分别由 no-child 项和 finite-child
链项支付，不能期待仅靠任一机制完成一般下界。

## 一个经基本词检查的混合候选

对任一词，令 \(\mathcal K(\sigma)\) 是所有彼此顶点不交 clean-chain families 的
集合，并定义严格可计算的最大链证书
\[
 C_{\rm chain}(\sigma):=
 \sup_{\mathcal F\in\mathcal K(\sigma)}
 \sum_{\mathcal C\in\mathcal F}\frac{2M_{\mathcal C}^3}{9h_{\mathcal C}N}.
\tag{6}
\]
R055 的局部引理立即给出 \(C_{\rm chain}(\sigma)\le\Psi_{\rm fin}(\sigma)\)。

下列命题是新的**待证候选**，不是本轮结论：
\[
 \boxed{\qquad
 \Psi_\infty(\sigma)+C_{\rm chain}(\sigma)\ge cN\log n
 \quad\hbox{for every one-peak word on }P_n\square P_n.
 \qquad}
\tag{MC}
\]
它通过两个最低限度的基本检验：

* 逐行词由 (3) 单独满足量级；
* comb 词由既有的 clean-chain 调度证书支付更强的 \(\Omega(n^3)\) 量级。

R058 只否定去掉 \(\Psi_\infty\) 后的纯链版本，不能否定 (MC)。

## 第一处未证步骤

当前没有从方格连通增长序的二维结构推出 (MC)。尤其不能简单按
\(b_p=N+1\) 与 \(b_p\le N\) 的点数分类：

* 无 sole-child 点可以在很早的 \(\ell_p\) 集中出现，此时
  \((l_p-1)/(N-l_p)\) 只有小量；
* 有 finite child 的点可以形成大量短、分叉或重叠链，R055 的 disjoint chain
  证书不一定捕获其总贡献；
* 将这些链重叠使用需要一个新的全局重数预算。R058 的逐行词说明不存在无预算的
  clean-chain 覆盖。

因此 (MC) 只是比 R055 更贴合两种已知基本机制的候选。一个真正推进的下一步是：
为任何 no-child 低贡献区域构造有限-child 的链证书，或证明这种区域的 rank profile
本身给出 (3) 型调和下界；目前两种替代之间没有已证的二维 bridge。

## 有限审计范围

<code>audit_mixed_terms.py</code> 以整数/有理数独立分解一个既有 \(8\times8\)
低 \(\Psi\) 搜索词、\(9\times9\) 逐行词和 \(5,7,9\) 阶 comb 词的 (1)，并对各
comb 核对 (5) 的精确有理上界。
它仅核对定义与这些指定词，不为 (MC) 提供一般证据。
