# R033：远峰纤维的正矩近似与端点缺口

## 设定

固定 $\varepsilon>0$，令 $E_{n,\varepsilon}$ 为两个峰的原图曼哈顿距离
$d>\varepsilon n$ 的事件。把 uniform top-two 切片按首次合并纤维 $f$ 分解；
$B_f$ 是该纤维的 top-two 词数，$q_f=a_f+b_f$ 是合并前两分支大小之和。
已有几何事实给出

$$q_f\ge d_f.\tag{1}$$

以下先限于非 singleton 纤维 $a_f,b_f\ge2$。写

$$
X_f=a_f-2,\qquad Y_f=b_f-2,\qquad M_f=X_f+Y_f=q_f-4,\qquad
\theta_f=X_f/M_f.
$$

这里 $\theta_f$ 只在 $M_f\ge1$ 时定义。对固定 $j\ge1$，仅在
$M_f\ge j$ 时使用下式。固定纤维的已证公式是

$$
p_{f,j+2}=p_{f,2}\mu_{f,j},\qquad
\mu_{f,j}=\frac{(X_f)_j+(Y_f)_j}{2(M_f)_j},\qquad
p_{f,2}=\frac{(X_f+1)(Y_f+1)}{(M_f+2)(M_f+1)}.\tag{2}
$$

这里 $p_{f,t}$ 是该纤维在 top-two 词中 $L=t$ 的比例，$L$ 为项目既定的
首邻点间隙。

## 不放回—放回耦合

令 $j$ 个带颜色球从含 $X$ 个红球、$Y=M-X$ 个蓝球的总体中不放回抽取；
$\mu_j$ 正是“全红或全蓝”的概率的一半和。用放回抽取时，相应量是

$$r_j(\theta)=\frac{\theta^j+(1-\theta)^j}{2}.\tag{3}$$

将有放回样本在没有重复时与不放回样本相同，并在出现第一次重复时重新抽取
不放回样本。此耦合的失配概率至多有放回样本中一对位置重合的概率，故并集界给

$$
\boxed{\quad |\mu_j-r_j(\theta)|\le\frac{j(j-1)}{2M},\qquad
0\le X\le M,\ 1\le j\le M.\quad}\tag{4}
$$

常数不依赖 $X$，包括 $X=0,M$ 的端点。

令 $\mathcal F_{n,\varepsilon}^{\rm g}$ 是远峰且 $a,b\ge2$ 的纤维，定义正测度

$$
\nu_{n,\varepsilon}=\sum_{f\in\mathcal F_{n,\varepsilon}^{\rm g}}
B_fp_{f,2}\,\delta_{\theta_f},\qquad
\eta_{n,\varepsilon}=
\frac12\sum_fB_fp_{f,2}(\delta_{\theta_f}+\delta_{1-\theta_f}).\tag{5}
$$

因此 $\eta$ 是 $[0,1]$ 上的正对称有限测度，且它的第 $j$ 阶 Hausdorff 矩为

$$m_j(\eta)=\int_0^1t^j\,d\eta(t)=\sum_fB_fp_{f,2}r_j(\theta_f).$$

若 $h_{t,E}^{\rm g}$ 表示这些非 singleton 远峰纤维对短桶 $h_t$ 的贡献，则 (2)--(4)
给

$$
\left|h_{j+2,E}^{\rm g}-m_j(\eta_{n,\varepsilon})\right|
\le\frac{j(j-1)}{2M_*}\,\nu_{n,\varepsilon}([0,1])
=\frac{j(j-1)}{2M_*}h_{2,E}^{\rm g},\tag{6}
$$

其中，只要 $M_*:=\lfloor\varepsilon n\rfloor-3\ge\max\{1,j\}$，由 (1) 可取

$$M_f\ge d_f-4\ge\lfloor\varepsilon n\rfloor-3=M_*.$$

所以每个固定 $j\ge1$ 下，远峰非 singleton 短桶确实渐近于同一个**正**对称
Hausdorff 矩列；误差为 $O_j((\varepsilon n)^{-1}h_{2,E}^{\rm g})$。有限 $M$ 的
$\mu_{f,j}$ 本身不在这里被当作正矩列。$j=0$ 单独是恒等式
$h_{2,E}^{\rm g}=\nu_{n,\varepsilon}([0,1])$，不需要定义 $\theta_f$。

### Singleton 分支

若 $\min(a_f,b_f)=1$，根方向平均后 $L$ 只在 $1$ 和 $q_f-1$ 取值，各占一半。
因此 $q_f=j+3$ 时它恰对 $h_{j+2}$ 贡献一半；只有 $q_f\ge j+4$ 才对该桶为零。
对固定 $j$，若 $\lfloor\varepsilon n\rfloor+1\ge j+4$，远峰的 (1) 保证后一条件，
它才不进入 (5) 中对应的短桶；并且 $p_{f,2}=0$。这不是可忽略的误差：
它可能对 $B_E$ 和尤其 $A_{2,E}$ 有大量贡献，见下节。

## 真正等价的正积分表示

已证稳定提升恒等式在每个几何事件上给

$$A_{2,E}=\sum_{t\ge1}t\,h_{t,E},\qquad B_E=\sum_{t\ge1}h_{t,E}.\tag{7}$$

对一个非 singleton 纤维，令 $m=M+2$、$x=X+1$、$y=Y+1$。随机交错中低根首邻点
的条件均值为

$$
e_f=\mathbb E_fL=\frac{m+1}{2}\left(\frac1{x+1}+\frac1{y+1}\right)
=\frac{M+3}{2}\left(\frac1{X+2}+\frac1{Y+2}\right).\tag{8}
$$

把 (5) 视为 $(M,X)$ 上的正离散测度
$\Lambda_E=\sum_{f\in\mathcal F_E^{\rm g}}B_fp_{f,2}\delta_{(M_f,X_f)}$，并置

$$
\alpha(M,X)=\frac{(M+2)(M+1)}{(X+1)(M-X+1)},\qquad
\beta(M,X)=\alpha(M,X)\frac{M+3}{2}
\left(\frac1{X+2}+\frac1{M-X+2}\right).\tag{9}
$$

则以下不是近似而是精确正积分：

$$
B_E^{\rm g}=\int\alpha\,d\Lambda_E,\qquad
A_{2,E}^{\rm g}=\int\beta\,d\Lambda_E.\tag{10}
$$

若 $\mathcal F_E^{\rm s}$ 是 singleton 纤维，则同样精确地

$$
B_E^{\rm s}=\sum_{f\in\mathcal F_E^{\rm s}}B_f,\qquad
A_{2,E}^{\rm s}=\frac12\sum_{f\in\mathcal F_E^{\rm s}}q_fB_f.\tag{11}
$$

于是

$$
\text{若 }A_{2,E}>0,\qquad \frac{B_E}{A_{2,E}}=
\frac{B_E^{\rm s}+\int\alpha\,d\Lambda_E}
{\frac12\sum_{f\in\mathcal F_E^{\rm s}}q_fB_f+\int\beta\,d\Lambda_E}.\tag{12}
$$

这给出所求的真正等价正测度判据：当 $A_{2,E}>0$ 时，
$B_E/A_{2,E}\ge c$ 当且仅当 (12) 的分子至少为 $c$ 倍分母。若 $E$ 为空，
只使用 (10)--(11) 的未归一化等式。它们不丢失 singleton 或端点质量。

若 $X/M\to\theta\in(0,1)$，则

$$
\alpha(M,X)\to\frac1{\theta(1-\theta)},\qquad
\beta(M,X)\to\frac1{2\theta^2(1-\theta)^2}.\tag{13}
$$

因而若想**仅由** (5) 的 $\theta$-矩极限恢复 (12)，一组充分条件是：singleton 项在
相关归一化下可忽略，且 $\alpha,\beta$ 对 $\Lambda_E$ 在 $\theta=0,1$ 附近一致可积。
这不是 (6) 的推论，也不是处理原题的必要前提；别的路线可以显式保留非零 singleton
或端点缺陷。

一个具体警告是 $X=0,M\to\infty$：所有固定 $j\ge1$ 的极限矩为 $1/2$，但
$p_2=1/(M+2)$、$\alpha=M+2$、$\beta\asymp M^2/4$。故固定阶正矩信息既不能恢复
$B_E$，更不能恢复 $A_{2,E}$。singleton 纤维更极端：它在所有固定的短桶
$h_{j+2}$ 中消失，却有 $\mathbb E L=q/2$。

## 与原题的准确关系

对远峰事件 $E=E_{n,\varepsilon}$，原条件律恰为

$$\mathbb P_2(E)=A_{2,E}/A_2.\tag{14}$$

(6) 给出的是远峰 top-two 短桶的渐近正矩表示；当 $B_E>0$ 时，(12) 给出远峰分子
内部 $B_E/A_{2,E}=1/\mathbb E_{B_E}L$ 的等价正积分表示。若 $B_E=0$，则
$A_{2,E}=0$。要从 top-two 远峰占比推到 (14)，在 $B_E>0$ 时有

$$\mathbb P_2(E)=\frac{B_E}{B}\,
\frac{\mathbb E_{B_E}L}{\mathbb E_BL}.\tag{15}$$

因此固定阶短桶关系、甚至其全部逐阶极限，不能单独证明原题的远峰概率趋零；
对于“只取 $\theta$ 矩极限”的路线，第一处未证步骤是上述 singleton/端点控制；
也可改用显式缺陷项或其他全局尾控制。

## 有理数核验

`check_mu_coupling.py` 用 `Fraction` 精确核验若干 $M,X,j$，包括 $X=0$、端点附近
与平衡情形；所有案例均满足 (4)。这只核验常数和实现，不代替耦合证明。
