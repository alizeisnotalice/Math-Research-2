# 任意二次平移保持外质量势

令 \(S\) 是布尔种子 \(G\) 的真实 Walsh 谱，\(e=(1,\beta)\) 的张成为 \(V\)。对每个真子空间 \(K<V\) 记
\[
m_G(K)=\sum_{e\notin K}|\widehat G(\beta)|,\qquad
L(G)=\max_{K<V}\frac{\operatorname{codim}_V K+1}{m_G(K)}.
\]
给定任意新地址维数，令 \(F(u,x)=(-1)^{c(u)}G(x+T(u))\)，其中 \(c,T_i\) 至多二次。则
\[
\boxed{L(F)\le L(G).}
\]
不要求旧支撑满谱、系数等模、极化阵共线或各频带极化秩相同。后续阶段可以让平移作用于任何旧数据和旧地址坐标。证明实际只用归一化权重 \(\nu_e=L(G)|\widehat G(\beta)|\) 所满足的割条件
\[
 \sum_{e\notin K}\nu_e\ge\operatorname{codim}_V K+1\qquad(K<V).
\tag{1}
\]

## 线性参数的矩形矩阵族

设 \(e\mapsto B_e\) 是从 \(V\) 到 \(\operatorname{Hom}(Y,X)\) 的线性映射，实际 \(e\) 张满 \(V\)。假定所有 \(B_e\) 的联合像张满 \(X\)，且 \(\bigcap_e\ker B_e=0\)。令 \(a=\dim X,\ b=\dim Y,\ r_e=\operatorname{rank}B_e\)，其中 \(a,b>0\)。则
\[
 \sum_e\nu_e(2^{r_e}-1)\ge a+b,\qquad
 \sum_e\nu_e(2^{r_e}-2^{-r_e})\ge a+b+1.
\tag{2}
\]
两函数 \(f(r)=2^r-1\)、\(g(r)=2^r-2^{-r}\) 都递增且超可加。对 \(g\) 可令 \(x=2^r,y=2^s\)，其超可加差为 \((x-1)(y-1)(1-1/(xy))\)。

按 \(a+b\) 归纳。若存在 \(y\ne0\)，使 \(X_0=\operatorname{span}_e B_e y\) 是 \(X\) 的真子空间，则 \(X_0\ne0\)。设
\[
 Y_0=\{z\in Y:B_e z\in X_0\text{ 对所有 }e\}.
\]
它含 \(y\)，并且严格小于 \(Y\)，因为原联合像为整个 \(X\)。限制族 \(B_e|_{Y_0}:Y_0\to X_0\) 仍有满联合像（仅 \(y\) 的像已经生成 \(X_0\)）且无公共核。商族 \(\overline B_e:Y/Y_0\to X/X_0\) 同样有满联合像且无公共核。对每个 \(e\)，
\[
 \operatorname{rank}B_e\ge\operatorname{rank}(B_e|_{Y_0})+
                                  \operatorname{rank}\overline B_e ,
\]
因为限制的像在 \(X_0\) 内。故可对两个更小的矩形族应用归纳，再用超可加性支付两种费用。第二种费用还得到两个单位的余额。

否则每个 \(y\ne0\) 都有 \(\operatorname{span}_e B_e y=X\)。在使用此情形之前可先将全体矩阵转置，使 \(a\ge b\)。若转置后的矩阵族有前述真子空间分割，则使用前一段。转置保持两个联合维数与各秩。由 (1) 对求值映射 \(e\mapsto B_e y\) 的核给出
\[
 \sum_e\nu_e\mathbf1_{\{B_e y\ne0\}}\ge a+1\qquad(y\ne0).
\]
在 \(Y\) 中均匀平均，得到
\[
 \sum_e\nu_e(1-2^{-r_e})\ge(a+1)(1-2^{-b}).
\]
对正秩 \(r\)，有 \(f(r)=2^r(1-2^{-r})\ge2(1-2^{-r})\)，且 \(g(r)=(2^r+1)(1-2^{-r})\ge3(1-2^{-r})\)；零秩的两式均为等式。于是 (2) 由
\[
 2(a+1)(1-2^{-b})\ge a+b,\qquad
 3(a+1)(1-2^{-b})\ge a+b+1
\]
得出。两式对固定 \(b\ge1\) 只须核对 \(a=b\)，因为差分别随 \(a\) 非递减。第一式归结为 \(2^b\ge b+1\)。第二式在 \(b=1\) 取等，在 \(b\ge2\) 由 \(2^b\ge4\) 直接得出。这完成归纳。

## 任意交替矩阵族

设 \(e\mapsto A_e\) 是由 \(V\) 线性参数化的交替矩阵族，联合像张满 \(D\)，其中 \(\dim D=t_0>0\)。记 \(\rho_e=\operatorname{rank}A_e\)。则
\[
 \sum_e\nu_e(2^{\rho_e/2}-1)\ge t_0,\qquad
 \sum_e\nu_e(2^{\rho_e/2}-2^{-\rho_e/2})\ge t_0+1 .
\tag{3}
\]
这里对 \(t_0\) 归纳。任选一个非零 \(x\in D^*\)，令 \(S_x=\operatorname{span}_e A_e x\)。因联合像张满且 \(A_e\) 对称，\(S_x\ne0\)。又交替性给 \(S_x\subset\ker x\)。取 \(z\in D\) 使 \(x(z)=1\)，再在 \(\ker x\) 内取 \(S_x\) 的余空间 \(H\)。把 \(A_e\) 的两侧投到 \(H\)，记其矩阵为 \(A_e^+\)。令 \(D_+\) 为这些投影矩阵的联合像，并将 \(H=D_+\oplus N\)。

在 \(D=D_+\oplus(N+\langle z\rangle)\oplus S_x\) 中，原矩阵为
\[
 A_e=\begin{pmatrix}
 A_e^+&0&C_e\\
 0&0&B_e\\
 C_e^{\mathsf T}&B_e^{\mathsf T}&D_e
 \end{pmatrix}.
\]
在 \(N\) 上的行没有 \(H\) 列，因为投影联合像位于 \(D_+\)。在 \(z\) 上的行同样只含 \(S_x\) 列，因为其原坐标泛函就是 \(x\)，而 \(A_ex\in S_x\)。这证实分块中的零。矩形块 \(B_e:S_x^*\to N+\langle z\rangle\) 的联合像张满 \(N+\langle z\rangle\)，因为完整原像张满 \(D\)，投影到这一坐标块只来自 \(B_e\)。其联合行空间也张满 \(S_x\)，因为各 \(z\) 行恰给 \(A_ex\)，这些行已张满 \(S_x\)。

对任意 \(e\)，写 \(k_e=\operatorname{rank}B_e\)。在矩形块中选 \(k_e\) 对配对方向，使对应交叉块为单位阵。由主子块及 Schur 补，
\[
 \rho_e\ge\operatorname{rank}A_e^++2k_e.
\tag{4}
\]
具体地，所选 \(2k_e\) 阶配对块形如 \(\left(\begin{smallmatrix}0&I\\ I&E\end{smallmatrix}\right)\)，其逆的右下角为零。由于跨到 \(D_+\) 的系数左半边为零，Schur 补仍为 \(A_e^+\)，所以两秩相加。

若 \(D_+\ne0\)，对它应用维数归纳。对 \(B_e\) 应用 (2)。若 \(D_+=0\)，只用 (2)。由 (4) 及两费用对半秩的超可加性，第一种费用支付至少 \(\dim D_++\dim(N+\langle z\rangle)+\dim S_x=\dim D\)。第二种至少再支付一个单位。这里 \(D_+\) 总是严格降维，因为 \(\dim S_x\ge1\) 且还去了一维 \(z\)。由此得 (3)。

## 返回父谱的任意切割

固定旧频率 \(\beta\) 后，父带系数是 \(\widehat G(\beta)\) 乘以 \((-1)^{c+\beta\cdot T}\) 的地址 Walsh 系数。令 \(M_e\) 为此相位的交替极化阵，则 \(M_e\) 随齐次列 \(e\) 线性。对二次相位 \(q\)，其符号函数沿 \(h\) 的自相关仅在 \(\ker M_e\) 上非零，值为 \((-1)^{q(h)+q(0)}\)。再作 Walsh 变换表明系数平方在与该线性符号相容的一个 \(\operatorname{im}M_e\) 陪集上恒为 \(2^{-\operatorname{rank}M_e}\)，其余为零。在公共根 \(\bigcap_e\ker M_e\) 上，\(q(h)+q(0)\) 关于 \(h\) 线性，且随 \(e\) 线性。把这些线性型延拓到整个地址空间，可用统一的齐次频率剪切将所有频带放在 \(W=(\bigcap_e\ker M_e)^\perp=\sum_e\operatorname{im}M_e\) 上。带内差生成全部 \(W\)，旧投影张满 \(V\)，故父齐次谱张成恰为 \(V\oplus W\)。频率剪切保留 \(L\)。

设父端真子空间为 \(X\)。写 \(U=\pi_VX,\ J=X\cap W,\ Q=W/J,\ u=\operatorname{codim}_VU\)。取 \(\pi:W\to Q\)，定义 \(A_e=\pi M_e\pi^{\mathsf T}\)，并令 \(Q_A\) 为其联合像。选择 \(Q=Q_A\oplus N\)。记 \(\rho_e=\operatorname{rank}A_e\)，\(k_e\) 为单侧投影 \(\pi M_e\) 再投到 \(N\) 的秩。在 \(W=Q_A\oplus N\oplus J\) 中，\(N\) 行只含 \(J\) 列。对这一交叉块使用与 (4) 相同的 Schur 补。另一方面，投影像在 \(Q_A\) 内包含 \(\operatorname{im}A_e\)，而在模 \(Q_A\) 后有 \(k_e\) 维。得到
\[
 \operatorname{rank}M_e\ge \rho_e+2k_e,\qquad
 p_e:=\operatorname{rank}(\pi M_e)\ge \rho_e+k_e .
\]
每个父外带的归一化费用 \(g_e\) 至少如下：
\[
\begin{aligned}
e\notin U:\quad
g_e&=2^{\operatorname{rank}M_e/2}
 \ge 1+(2^{\rho_e/2}-1)+(2^{k_e}-1),\\
e\in U:\quad
g_e&\ge2^{\operatorname{rank}M_e/2}(1-2^{-p_e})
 \ge (2^{\rho_e/2}-2^{-\rho_e/2})+(2^{k_e}-1).
\end{aligned}
\]
用 \(k_e\) 作费用，令 \(B\) 为按费用递增选出的一个贪心基。阈值 \(\lambda\ge0\) 下，费用不超过 \(\lambda\) 的列所张成的空间，其余维恰为 \(B\) 中费用超过 \(\lambda\) 的列数。由 (1) 得大于 \(\lambda\) 的总质量至少为此余维加一，只在余维非零时计入这一单位。对阈值积分便得 \(\sum_e\nu_e k_e\ge\sum_{e\in B}k_e+\max_{e\in B}k_e\)。由于任一原基的单侧矩阵像都张满 \(N\)，在 \(N\ne0\) 时，\(\sum_e\nu_e(2^{k_e}-1)\ge\dim N+1\)。对两侧投影 \(A_e\) 使用 (3)。

若 \(U<V\)，外带基线支付 \(u+1\)，两侧及单侧支付至少 \(\dim Q_A+\dim N\)，因为内带的 \(\phi\) 费用不小于 \(h\) 费用。若 \(U=V\)，所有标签都使用第二行。\(N\ne0\) 时单侧已含一个额外单位。\(N=0\) 时，因 \(X\) 为真子空间，必有 \(Q_A=Q\ne0\)，由 (3) 的第二式支付这个单位。综上
\[
 L(G)\,m_F(X)\ \ge\ \sum_e\nu_e g_e
 \ \ge\ u+\dim Q+1=\operatorname{codim}_{V\oplus W}X+1.
\]
故 \(L(F)\le L(G)\)。

特别地，从任意固定种子出发，允许任意深度、任意新地址维数、每代任意混合极化秩的这类递归均有 \(P(F)\le L(G)\)。从字符出发上界为二。这不涵盖第85轮中把已经存在的控制位重用于同一层的三角门，此处每一次平移使用独立新地址。一般 A28 的外质量势仍无绝对上界。

[第82轮逐带公式](082-phases.md) · [第85轮单侧投影的 Schur 补](085-schur.md)
