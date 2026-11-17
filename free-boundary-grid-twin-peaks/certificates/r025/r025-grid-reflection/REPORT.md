# r025：方格双峰根的反射—协方差路线审计

日期：2026-11-16。对象始终是自由边界方格
\(G_n=P_n\square P_n\)，在均匀标号的条件律 \(K_n=2\) 下把无序峰对
均匀定向为 \((X,Y)\)。本目录是临时研究记录，不修改项目仓库。

## 结论

截至本轮，**没有得到**全 \(n\) 的宏观负协方差，也没有得到
\(\operatorname{Var}(X_x)\ge c n^2\)。但首次合并表示把两件事各自化成了一个
完全精确、带完成数权重的分支矩不等式。整条竖直反射保持该协方差求和的每一项
不变，因而单靠全局反射对称不能给出符号；反射若要起作用，必须是一个受控的
**单支**交换，同时比较其完成数权重。该比较正是目前第一处未证步骤。

已有严格有限反例还排除了三种较强、但不必要的局部替代品：逐根镜像迁移、固定
paired-label 反射轨道上的 split 偏向，以及所有竖切交叉协方差逐项非正。
这不否定全局负协方差，更不构成原方格问题的反例。

## 1. 首次合并下的精确坐标矩

固定横坐标，并记

\[
t(v)=2x(v)-(n-1).
\]

令 \(F=(S,w)\) 跑过首次合并事件，\(A_F,B_F\) 是加入 \(w\) 前的两个
分支，\(c_F=h_G(S)\binom{|S|-1}{|A_F|}\) 是已有的正完成数权重。为避免把
分支 \(A_F\) 和总计数混淆，定义

\[
 \alpha(H)=\sum_{u\in H}a_H(u),\qquad
 m(H)=\sum_{u\in H}a_H(u)t(u),\qquad
 v(H)=\sum_{u\in H}a_H(u)t(u)^2.
\]

其中 \(a_H(u)\) 是诱导图 \(H\) 中以 \(u\) 为唯一峰的相对标号数。将
`weighted-root-covariance.md` 的精确根核逐项求和，得到

\[
 A_2(n)=\sum_Fc_F\alpha(A_F)\alpha(B_F),\tag{1}
\]
\[
 \mathbb E_2[t(X)t(Y)]
   ={1\over A_2(n)}\sum_Fc_Fm(A_F)m(B_F),\tag{2}
\]
\[
 \mathbb E_2[t(X)^2]
   ={1\over2A_2(n)}\sum_Fc_F
   \{v(A_F)\alpha(B_F)+\alpha(A_F)v(B_F)\}.\tag{3}
\]

证明只须展开根核：在固定 \(F\) 时，两个根分别落在两个分支，权重为
\(a_{A_F}(p)a_{B_F}(q)\)；对 \(p,q\) 求和即分别产生 \(m\) 或 \(v,\alpha\)。
随机定向的两个顺序对 (2) 的乘积相同，对 (3) 则给出因子 \(1/2\)。这也核对了
归一化。

D4 对称给出 \(\mathbb E_2t(X)=0\)，所以

\[
 \operatorname{Cov}_2\!\left(X_x-{n-1\over2},Y_x-{n-1\over2}\right)
 ={1\over4A_2(n)}\sum_Fc_Fm(A_F)m(B_F),\tag{4}
\]
\[
 \operatorname{Var}_2\!\left(X_x-{n-1\over2}\right)
 ={1\over8A_2(n)}\sum_Fc_F
 \{v(A_F)\alpha(B_F)+\alpha(A_F)v(B_F)\}.\tag{5}
\]

因此，负协方差的**精确等价候选**是

\[
 \boxed{\ \sum_Fc_Fm(A_F)m(B_F)\le0.\ }\tag{C}
\]

而横坐标方差 \(\ge\kappa n^2\) 的精确充分且必要条件是 (5) 右端分子
\(\ge8\kappa n^2A_2(n)\)。它们都是对原条件律的陈述，不引入新的模型。

## 2. 为什么整条反射本身不给 (C) 的符号

设 \(R\) 是完整方格的竖直反射。它把 \(F\) 双射到 \(RF\)；分支可能互换，
但二项式权重及完成数不变。且

\[
 \alpha(RH)=\alpha(H),\qquad v(RH)=v(H),\qquad m(RH)=-m(H).
\]

故 (C) 中 \(F\) 的项在 \(RF\) 下变为

\[
c_{RF}m(RA_F)m(RB_F)=c_Fm(A_F)m(B_F),
\]

而不是其相反数。也就是说，反射轨道配对只是把同一符号项重复一次；它既不能
消去正项，也不能压过正项。这是一个严格的结构性限制，而非有限计算现象。

若只反射一个分支，则一般既不保持首次合并状态：\(RA_F\) 可能和 \(B_F\) 相交，
也未必仍仅由同一个 \(w\) 首次连接。即使偶然仍合法，也没有理由使
\(h_G(S)\)（从而 \(c_F\)）保持不变。故一个可行的 cluster-exchange 证明首先
必须建立一个带权注入/运输，足以比较 (C) 的正、负两部分；不能只证明几何形状
被镜像。

这给出本轮的第一未证条件：构造或估计一个**单支反射的累计完成数比较**，使
\[
 \sum_{F:m(A_F)m(B_F)>0}c_Fm(A_F)m(B_F)
 \le
 \sum_{F:m(A_F)m(B_F)<0}c_F|m(A_F)m(B_F)|.
\tag{6}
\]
式 (6) 与 (C) 等价，因而不是新假设；它明确指出交换必须同时处理的权重和方向。
目前没有此比较。

## 3. 已核验的障碍及其强度

- 固定一个根再把另一个根镜像向内移动的点态比较已在 \(3\times3\) 失败：
  \(K((0,0),(2,2))=7872>2688=K((0,0),(0,2))\)，尽管距离从 4 降为 2。
- 在 \(4\times4\)，固定 paired-label 反射 orbit 的 256 个朝向中，恰双峰的
  同侧/异侧数可为 12/6。因此不能逐 orbit 注入 same 到 split。
- 精确切协方差在 \(6\times6\) 有
  \(\operatorname{Cov}(s_0(X),s_4(Y))>0\)；所以不能把坐标协方差写成切矩阵后
  再逐项证明非正。
- 中央切的有效低峰—镜像点交换在有效域上是单射，但其必要条件已经包含额外的
  标签不等式，\(4\times4\) 有同侧双峰标号不在有效域。该法当前仅给
  \(\mathbb P(\mathrm{split})\ge\mathbb P(E)\)，不能控制 \(E\) 的条件质量，
  也不能给 (C) 的坐标权重。

这些结论的来源和哈希列在下面的记录中。它们排除的是特定加强命题；全局的 (C)
仍可能为真。

## 4. 一峰转移可怎样服务于方差，而尚缺什么

令原双峰根的未归一化边际为

\[
 r_n(p)=\sum_{q\ne p}K_n(p,q)
 =\sum_Fc_F\bigl[
 {\bf1}_{p\in A_F}a_{A_F}(p)\alpha(B_F)+
 {\bf1}_{p\in B_F}a_{B_F}(p)\alpha(A_F)\bigr].\tag{7}
\]

于是 \(\mathbb P_2(X=p)=r_n(p)/(2A_2(n))\)。对任意固定
\(0<\delta<1/2\)，若能证明中央竖带的根边际反集中

\[
 {1\over2A_2(n)}\sum_{|x(p)-(n-1)/2|<\delta n}r_n(p)\le\eta<1\tag{8}
\]

（常数对充分大的 \(n\) 一致），便立即有
\[
 \operatorname{Var}_2(X_x)\ge\delta^2(1-\eta)n^2.\tag{9}
\]

这是一条干净的“根转移 \(\Rightarrow\) 宏观方差”目标。不过已知
\(a_x(H)=\sum_{y\sim x}h_H(\{x,y\})\) 的可逆一峰链只控制**固定**诱导图 \(H\)
的根；(7) 混合了所有随机分支 \(A_F,B_F\)，并加上完成数 \(c_F\)。把固定图的
局部 conductance 下界直接加总为 (8) 的加权不等式，正是另一处尚未证明的步骤。
已有删除递推还必须处理带洞/细长的任意连通诱导子图，故不能只对完整方格归纳。

因此一峰根转移目前提供了正确的接口，却没有证明方差下界。仅有 D4 对称也允许
根边际集中在中央 \(O(1)\) 宽区域，故不可能单独推出 (9)。

## 5. 本轮实际有限核验

`audit_n3_root_moments.py` 完整枚举 \(9!\) 个 \(3\times3\) 标号，独立核对上述
有序根矩的因子：

\[
A_2(3)=143808,\quad
\mathbb E[t(X)t(Y)]=-170/107,\quad
\mathbb E[t(X)^2]=2110/749.
\]

前者给单坐标协方差 \(-85/214\)，与已有根核表相符；两个反射核数也复得
7872 和 2688。它只检查有限归一化，不是全 \(n\) 证据。

按用户指定的远程计算流程，先尝试 `ssh macbook-math` 作轻量环境核验；在 5 秒
内 banner 超时，因而没有在远端启动或遗留进程。上述小穷举在本机 Python 3 完成，
耗时约 1.1 秒；检查后无本任务计算进程存活。

## 6. 对原问题的含义

若未来既证明 (C)，又从 (5) 或 (8) 得到
\(\operatorname{Var}(X_x)\ge\kappa n^2\)，则 D4 给出
\(\mathbb E\|X-Y\|_2^2\ge4\kappa n^2\)，从而提供远峰正质量的路线。反过来，
本轮没有证明 (C)、(8) 或任何宏观下界，不能宣布原方格命题被证伪或被证明。

## 输入、复现和来源记录

本轮新文件：

```sh
python3 /Volumes/data/CodexTemp/T2/a022-grid-reflection-r025/audit_n3_root_moments.py
```

输出为 `audit_n3_root_moments.json`。读取的既有资料（SHA-256）为：

```
986a2e1f0d806be9ba86fe885fa4d8f9a57b9878fc15e82fd457d7a16755d6f8  proofs/weighted-root-covariance.md
0892b39043d639543123b5d4f896c545f6183dd61b9e90624156a0bcd1555f56  proofs/centered-cut-covariance.md
3487799415ad825c2c4d69f2c1b1a03fc45dd7758246b18e580fa982e35f07bb  proofs/central-cut-swap-audit.md
55411cd48d538112c9711f846c2fba1e740916ffebaca8dd204df71fdeaba86b  certificates/r020-reflection-pair/reflection-pair-orbit-audit.md
757828b7fd7cff404ec0cb72ae12142c3d4bf125b0e452859a6496c1d3f79612  certificates/r021-root-transport/REPORT.md
```

所有来源位于
`/Volumes/data/CodexProjects/T2/free-boundary-grid-twin-peaks/`；本轮未修改这些
来源、项目仓库或 Git 状态。
