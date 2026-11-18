# 仿射分割同带多原像的能量与符号预算

设 \(L=\bigcap_i L_i\) 为所有物理块共同方向，固定一个商标签 \(\beta\in L^*\)。在选定 \(X\cong (X/L)\oplus L\) 的线性坐标后，令 \(P=P_\beta\) 为该带在 \(X/L\) 上的系数函数。商仿射块 \(c_i+H_i\) 不交覆盖，质量为 \(p_i\)。令 \(P_i\) 为 \(P\) 在该块上的限制（按所选原点移到 \(H_i\)），\(m_i=s(P_i)\)，\(e_i=\|P_i\|_2^2\)。只在该带活动块上求和，记 \(E=\|P\|_2^2=\sum_i p_i e_i\)。

## 局部谱质量到全局支持

写 \(b_\xi=\widehat P(\xi)\)，\(c_{i,\eta}=\widehat P_i(\eta)\)，\(\pi_i\colon (X/L)^*\to H_i^*\) 为限制。反演后限制给精确恒等式
\[
c_{i,\eta}=\sum_{\xi:\pi_i(\xi)=\eta}b_\xi\,\chi_\xi(c_i).
\]
局部频率 \((i,\eta)\) 的完整提升是纤维 \(C_{i,\eta}=\pi_i^{-1}\{\eta\}\)。故真谱必命中每个活动纤维，且对非活动纤维不可能恰相交一个真谱点。任取非负纤维权重 \(w_{i,\eta}\)，若每个环境频率所处的加权活动纤维总权不超过1，则真谱大小至少为这些权重的总和。

同一恒等式对每块给 \(\|\widehat P_i\|_1\le\|b\|_1\le\sqrt{s_\beta E}\)。定义有效原像数
\[
\kappa_i=\frac{\|\widehat P_i\|_1^2}{e_i}\le m_i,\qquad
\rho_i=\frac{p_i e_i}{E},
\]
则
\[
s_\beta\ge\max_i\frac{\kappa_i e_i}{E}
\ge\sum_i\frac{\kappa_i\rho_i^2}{p_i}.
\]
局部该带系数平坦时 \(\kappa_i=m_i\)。若其活动块能量密度 \(e_i\) 相同，\(\rho_i=p_i/\sum_jp_j\)，即恢复第41轮想要的多原像频谱配重，而活动域无需仿射。若在不同标签上 \(e_{i,\beta}=a_i b_\beta\) 可按正数分离，同一全局指数次序改用速率 \(p_i a_i\) 即可取能量配重。不具这种分离性时，不能把不同标签的任意主块边际当作同一全局随机顺序。

另有物理幅值估计 \(s_\beta\ge E/(\sum_i p_i\|P_i\|_1)^2\)，不弱于仅用物理活动概率的不确定性。若已知逐频升谱强度，把 \(U_\xi=(\sum_i p_i|c_{i,\pi_i\xi}|)^2\) 从大到小排序，则前 \(s_\beta\) 项之和至少为 \(E\)。

## 只用活动局部谱符号的额外改进

取 \(\omega_i\ge0\) 且至少一个活动块权重为正。令 \(r_{i,\eta}=\operatorname{sgn}(c_{i,\eta})\)，零值的符号取零。由限制恒等式
\[
 \sum_i\omega_i\|\widehat P_i\|_1
 =\sum_\xi b_\xi \sum_i\omega_i\chi_\xi(c_i)r_{i,\pi_i\xi}.
\]
因此
\[
s_\beta\ge
\frac{(\sum_i\omega_i\|\widehat P_i\|_1)^2}
{E \left(\max_\xi|\sum_i\omega_i\chi_\xi(c_i)r_{i,\pi_i\xi}|\right)^2}.
\]
分母里的符号相消是几何位置与局部符号共同决定的。只取单块权重时退回该块的能量界。也可换成任意局部谱测试向量 \(u_{i,\eta}\)，把分子改成 \(|\sum u_{i,\eta}c_{i,\eta}|^2\)，分母改成相应线性组合最大值的平方。

## 确定性共享秩条件

沿用第41轮 \(N=\sum_i L_i^\perp\)、\(V_i=\operatorname{span}\{(\beta,1):\beta\in B_i\}\)。由于各块物理支撑不交，某个标签在局部出现则整带必非零，所以 \(\sum_i V_i\) 是全局谱增广张成的精确商投影。令 \(D_i=a_i+2\)，\(\gamma_i=D_i/s_i\)，并假定 \(D_i\le2\sqrt{s_i}\)。定义确定性秩剩余
\[
B_{\mathrm{rank}}
=\sum_iD_i-\left(\dim N+\dim\sum_i V_i+1\right).
\]
第2轮法向秩界和局部投影保证非负。任一共享带 \(\beta\)，取已证下界 \(\lambda_\beta\le s_\beta\)，允许来自上述能量、符号或其他估计。如果
\[
B_{\mathrm{rank}}\ \ge\
\sum_{\text{共享 }\beta}\left[
\sum_i\gamma_i m_{i\beta}
-\sqrt{\lambda_\beta\sum_i p_i\gamma_i^2m_{i\beta}}
\right],
\]
则 \(a+2\le2\sqrt{s}\)。

证明是在每个共享带分配 \(q_{i\beta}\ge0\)，满足 \(\sum_i q_{i\beta}^2/p_i\le\lambda_\beta\)。第 \(i\) 块对应标签的 \(m_{i\beta}\) 个局部原像格各填 \(q_{i\beta}/\sqrt{m_{i\beta}}\)，私有标签格仍填1。该向量的平方和 \(T_i\) 满足 \(\sum_i T_i/p_i\le s\)，因为私有带无跨块相消。局部向量的线性和为 \(s_i-\sum_{\rm shared}(m_{i\beta}-q_{i\beta}\sqrt{m_{i\beta}})\)。逐带最大化 \(\sum_i\gamma_i q_{i\beta}\sqrt{m_{i\beta}}\) 的值恰为上述根式，由 Cauchy--Schwarz 达成。于是秩剩余付清线性损失后，再用局部 Cauchy--Schwarz 和输入界得到 \(a+2\le2\sum_i\sqrt{T_i}\le2\sqrt{s}\)。

若不想计算全体 \(\sum_i V_i\)，可用确定性下界
\[
B_{\mathrm{rank}}\ge \sum_i (a_i+1-\dim V_i)+\dim\mathcal R
\]
其中 \(\mathcal R\le\bigoplus_i V_i\) 由各同标签增广向量在不同块的两坐标和张成。假定每块的共享增广标签线性独立，则每个共享标签选一个主块并把其余块连为星形关系，全部这些关系独立。在共享局部谱平坦，且每个 \(i,\beta\) 满足 \(\gamma_i m_{i\beta}\le1\) 时，已足够任意能量权重：因 \(\sum_i \gamma_i m_{i\beta}(1-\rho_{i\beta})\le\#\{i:\beta \text{ 活跃}\}-1\)。这种情形不需要同一全局随机顺序。

## 支持约束的精确线性判据

归一局部物理包的 Fourier 矩阵 \(A_{\xi,(i,\eta)}=\sqrt{p_i}\chi_\xi(c_i)1_{\xi\in C_{i,\eta}}\) 为一个实正交矩阵。若局部非零指数集为 \(I\)、候选全谱为 \(T\)，取 \(B=A^\top_{I^c,T}\)。要求存在某实带信号同时给出这两个支持图样，当且仅当每个 \(e_\xi\)（\(\xi\in T\)）与每个活动限制行 \(A^\top_{(i,\eta),T}\) 都不在 \(B\) 的行空间内。因为全谱系数应属于 \(\ker B\)，并可同时避开上述有限个真超平面。给定局部符号后可再加严格的实线性不等式；完整 Boolean 粘合还要求各数据带联合符合 Boolean 约束，单带线代没有代替这一点。


在多个原像的非仿射活动域上，必须按实际分割确定 $L$。[几何构造及完整局部预算](052-verified-fibers.md) 通过细分一个非活动地址的横向空间，使共同切向恰为所需的数据空间。

## 符号信息改变秩判据的无界族

取选择变量 $x\in\mathbb F_2^3$、$w\in\mathbb F_2^2$、$u\in\mathbb F_2^{2k}$，以及两个数据位 $y$。令 $Q(v)$ 为偶数个位依次两两配对后的乘积之和，令
$$h(x,w,u)=(-1)^{x_1x_2+x_1x_3+x_2x_3+Q(w)+Q(u)},\qquad
F=h\,(-1)^{y_1y_2}.$$
对每一对 $(u,w)$，取依 $w$ 从 $1,2,4,7$ 中选的三位向量 $n_w$，以 $n_w\cdot x=0,1$ 将该层分成两块。共有 $m=8\cdot4^k$ 个等质量块。因这四个法向张满三维，公共切向恰为 $y$ 空间。四个数据标签全部出现在每个块上。

三变量二次型的极化双线性型的核是 $\langle(1,1,1)\rangle$。四个 $n_w$ 对此核向量均取1，故该二次型限制到每个二维切片的极化形式均非退化。于是局部 $h$ 的四个 Walsh 系数等模 $1/2$，每个数据标签都在本块中有四个原像，局部总谱有16点且 $D_i=6$。全局三变量符号的谱有四点，绝对系数 $1/2$；这也可取核的一个补空间后，把非退化二变量型与核上的字符相乘直接得到。另两个偶维配对型给独立的平坦谱。故 $h$ 的真谱大小恰为 $16\cdot4^k$，模为 $1/(4\cdot2^k)$。

对一个固定的数据标签，数据 Fourier 系数为 $c=\pm1/2$，故该带能量 $E=c^2$。只用单块能量给 $\lambda=4$。在上述符号公式中取所有活动块权重为1。每个局部符号等于对应局部系数乘 $2/|c|$，全局合成公式遂给
$$\max_\xi\left|\sum_i \chi_\xi(c_i) r_{i,\pi_i\xi}\right|=2m\|\widehat h\|_\infty.$$
分子为 $(2m|c|)^2$，所以符号下界为 $1/\|\widehat h\|_\infty^2=16\cdot4^k$，恰好等于真带大小。此处 $c_i$ 表示物理块基点。

共同法向秩为 $5+2k$，共同标签增广秩为3，故秩扣减 $B_{\rm rank}=6m-(9+2k)$。四个标签具有相同下界 $\lambda$ 时，判据的右端为 $6m-3\sqrt\lambda$。单块能量给 $6m-6$，超过可用扣减；符号估计给 $6m-12\cdot2^k$，在每个 $k\ge0$ 均能支付。因此在同一秩判据中，局部符号改变了结论。该显式函数另有直接的乘积表示；此对比衡量指定分割上的证书强度。
