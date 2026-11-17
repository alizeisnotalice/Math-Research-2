# R061：前缀边界到 finite-child 父点的严格桥

## 定义

固定 \(P_n\square P_n\) 上的一峰词，\(N=n^2\)。令 \(U_t\) 是前 \(t\) 个点，
\[
 D_t=V\setminus\bigl(U_t\cup N(U_t)\bigr),\qquad
 B_t=e(U_t,V\setminus U_t).
\]
令 \(s_t\) 是 \(U_t\setminus\{r\}\) 中加入时恰有一个早邻点的顶点数。对 eligible
父点 \(p\)，仍记 \(b_p\) 为第一个以 \(p\) 为唯一早邻点的 child 的秩（无此 child
时 \(b_p=N+1\)），并记
\[
 F_t:=\{p\ {\rm eligible}:b_p\le t\},\qquad m_t=|F_t|.
\]

## 边界—sole-child 不等式（严格）

从根开始，加入顶点 \(v\) 时若有 \(e(v)\) 条早边，边界变化为
\(\deg(v)-2e(v)\)。当 \(e(v)=1\) 时它至多为 \(2\)，当 \(e(v)\ge2\) 时它至多为
零。因此
\[
 B_t\le\deg(r)+2s_t.
\tag{1}
\]
另一方面，每个已支配但未出现的顶点 \(q\in V\setminus(U_t\cup D_t)\) 至少附着
一条不同的边界边，故
\[
 B_t\ge N-t-|D_t|.
\tag{2}
\]
合并得
\[
 s_t\ge
\left\lceil\frac{N-t-|D_t|-\deg(r)}2\right\rceil.
\tag{3}
\]
若右端为负，按 \(s_t\ge0\) 理解。

每个父点至多拥有四个 singleton child，所以这些 \(s_t\) 个 child 至少涉及
\(\lceil s_t/4\rceil\) 个父点。根及其 \(\deg(r)\) 个邻点至多去掉
\(1+\deg(r)\) 个非eligible父点。故
\[
 \boxed{\quad
 m_t\ge
 \max\left\{0,\left\lceil\frac{s_t}{4}\right\rceil-(1+\deg(r))\right\}.
 \quad}
\tag{4}
\]
特别地，使用 \(\deg(r)\le4\) 的稍弱形式为
\[
 m_t\ge
\max\left\{0,\frac{N-t-|D_t|-4}{8}-5\right\}.
\tag{5}
\]

这正确排除以下可能性：在一个早已支配或近支配的线性大小前缀中，只有
\(o(N)\) 个 finite-child eligible 父点。

## 一个严格的加权推论

对 \(p\in F_t\)，有 \(R_p=b_p-l_p-1<t\)，故
\[
 \frac{l_p-1}{R_p}\ge\frac{l_p-1}{t}.
\]
此外，每个秩至多作为四个 \(l_p\) 出现，eligible 点有 \(l_p\ge2\)。若把
\(l_p-1\) 在 \(F_t\) 中排序，其第 \(j\) 个至少为 \(\lceil j/4\rceil\)。所以
\[
 \sum_{p\in F_t}\frac{l_p-1}{R_p}
\ge\frac1t\sum_{j=1}^{m_t}\left\lceil\frac j4\right\rceil
\ge\frac{m_t(m_t+4)}{8t}.
\tag{6}
\]
因 \(F_t\) 的每项皆属 \(\Psi_{\rm fin}\)，有
\[
 \boxed{\quad
 \Psi_{\rm fin}\ge\frac{m_t(m_t+4)}{8t}
 \quad(1\le t\le N).\quad}
\tag{7}
\]

例如，若对某个 \(t\le(1-\eta)N\)（\(\eta>0\) 固定）有
\(|D_t|=o(N)\)，则 (3)--(7) 至少强制 \(\Psi_{\rm fin}=\Omega_\eta(N)\)。
等价地，所需条件是 \(N-t-|D_t|\) 具有线性量级。这将 R061 的纯计数桥提升为
一个真实的 \(l_p/R_p\) 加权下界。

## 为什么这还不产生 \(N\log n\)

式 (7) 只能对一个 \(t\) 使用一次；集合 \(F_t\) 随 \(t\) 单调嵌套，不能把不同
时间的下界直接相加。更关键地，(3) 只控制 first singleton-child **数量**，而不控制
其 \(b_p-l_p\) 在不同尺度的分布。很多父点可以有很早的 \(l_p\)，或其 first child
分散在很晚的秩，从而使单一时间切片最多给线性量级。

因而当前第一缺口是一个带重数预算的时间分层命题：需要把
\(F_{2t}\setminus F_t\) 等互不相交的 first-child 秩层，与相应父点的 \(l_p\) 和
\(R_p\) 同时联系起来。边界不等式只给累计 \(s_t\)，不控制这些层的增量；
不允许把 (7) 对 dyadic \(t\) 相加。

## 有限审计

<code>audit_boundary_solechild.py</code> 对一个既有 \(8\times8\) 低 \(\Psi\) 词、
\(9\times9\) 逐行词和 \(9\times9\) comb 词，逐前缀以整数/有理数验证 (1)--(7)。
这是定义和常数的有限交叉核验；一般结论由上述边界增量和容量计数证明承担。
