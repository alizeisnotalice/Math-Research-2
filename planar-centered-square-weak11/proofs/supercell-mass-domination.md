# 超胞元质量支配输运与 (R=5) 尾界

日期：2026-11-15（Asia/Shanghai）。本稿是超胞元输运的质量支配版本。它为
“移动部分副本并调整其质量”的候选保留完整 halo 所需的全半径截断接口；不构造
候选、不计算其 halo，也不产生新的 (C_2) 下界。

## 1. 带标签的基线与目标超胞元

令
\[
 \mu=\sum_{i=1}^n\sum_{j,l\in\mathbb Z}
 w_i\delta_{p_i+(jP_x,lP_y)},\qquad
 p_i\in[0,P_x)\times[0,P_y),\quad w_i>0,
 \quad m=\sum_iw_i>0,
 \tag{1}
\]
并取 \(P_x,P_y,\lambda>0\)。选择整数 \(s_x,s_y\ge1\)，以
\(j=a+s_xk\)、\(l=b+s_yq\) 唯一分解原周期的平移指标，其中
\(0\le a<s_x\)、\(0\le b<s_y\)。记
\[
 z_{a,b,i,k,q}=p_i+((a+s_xk)P_x,(b+s_yq)P_y).
\]

令目标质量和位移满足
\[
 0\le w'_{a,b,i}\le\kappa w_i,\qquad
 \|h_{a,b,i}\|_\infty\le\delta,
 \tag{2}
\]
其中 \(\kappa>0\)、\(\delta\ge0\)，且并非所有 \(w'_{a,b,i}\) 都为零。定义
\[
 \nu=\sum_{a,b,i,k,q}w'_{a,b,i}
 \delta_{z_{a,b,i,k,q}+h_{a,b,i}}.
 \tag{3}
\]
它是 \(s_xP_x\mathbb Z\times s_yP_y\mathbb Z\)-周期测度。目标的实际单胞元质量是
\[
 m_\nu=\sum_{a=0}^{s_x-1}\sum_{b=0}^{s_y-1}\sum_iw'_{a,b,i}>0.
 \tag{4}
\]
它通常既不等于 \(s_xs_ym\)，也不等于 \(\kappa s_xs_ym\)。质量守恒是搜索时
可另加的约束，不是本引理的假设。

点坐标若跨越新胞元边界，程序可对 \([0,s_xP_x)\times[0,s_yP_y)\) 作模周期
约化并合并碰撞点；证明始终使用 (3) 的未约化标签坐标。故每个目标标签有唯一的
基线标签，wrap 和碰撞均不破坏 (2)。零目标质量只是不贡献相应标签；基线质量严格
为正使逐标签的比值/支配有定义。

## 2. 质量支配输运

对 \(Q(x,r)=x+[-r,r]^2\)，任意 \(x\in\mathbb R^2\)、\(r>0\) 有
\[
 \boxed{\quad\nu(Q(x,r))\le\kappa\mu(Q(x,r+\delta)).\quad} \tag{5}
\]

逐个标签求和即可证明。若 \(z+h\in Q(x,r)\)，则
\(z\in Q(x,r+\delta)\)；此目标标签贡献至多 \(\kappa w_i\)，正是右边相应基线
标签的 \(\kappa\) 倍。闭方形边界、共同半径、无穷平移索引及碰撞均由该逐标签包含
覆盖。局部有限性保证有界方形中只有有限个非零项。

式 (5) 不比较两边的极大函数：右边平均仍应除以目标半径的 \(4r^2\)，不能替为
\(4(r+\delta)^2\)。

## 3. 全尺度尾界与可认证的 \(\kappa\) 域

基线逐轨道计数给
\[
 \mu(Q(x,s))\le m\left(\frac{2s}{P_x}+1\right)
                       \left(\frac{2s}{P_y}+1\right).
 \tag{6}
\]
设 \(\rho=m/(P_xP_y)\)，并定义未乘质量支配常数的函数
\[
 G_\delta(r)=\rho\left(1+\frac{P_x+2\delta}{2r}\right)
                    \left(1+\frac{P_y+2\delta}{2r}\right).
 \tag{7}
\]
由 (5)–(6)，
\[
 \frac{\nu(Q(x,r))}{4r^2}\le \kappa G_\delta(r).
 \tag{8}
\]
正因子使 \(G_\delta\) 在 \(r>0\) 严格递减。因此固定 \(R>0\) 时，
\[
 \boxed{\kappa G_\delta(R)<\lambda}\tag{9}
\]
蕴含一切严格 \(\lambda\)-见证半径均小于 \(R\)。这里没有把算子人为限为
`r<R`。

固定 \(\delta\) 后，此尾界所允许的精确临界值是
\[
 \kappa_{\max}(\delta;R)=\frac{\lambda}{G_\delta(R)},\qquad
 0<\kappa<\kappa_{\max}(\delta;R).
 \tag{10}
\]
边界等号不通过此严格尾界验收；\(\kappa_{\max}\) 是可认证范围的上确界，不是
可直接接受的参数。该公式适用于任意 \(\lambda>0\)。

## 4. R8 基线的实际有理核验

冻结 R8 输入
`computation/inputs/r008_final_0611.json` 的 SHA-256 是
`6c2e5d2d8fba5063bd45ecd4467f8c74cca5c0edbef958f90146a103029d5910`。
它有 \(\lambda=1\)、\(R=5\)、
\[
 m=\frac{65941670400277227073796193886781311}
 {12500000000000000000000000000000000},\quad
 P_x=\frac{328734960629}{10^{11}},\quad
 P_y=\frac{175603458831370317}{64000000000000000}.
\]
对 \(\kappa=1001/1000\)，标准库 `Fraction` 的实际结果是：

| \(\delta\) | \(G_\delta(5)\) | \(\kappa G_\delta(5)<1\) |
| --- | --- | --- |
| \(1/1000\) | \(23828111711084880228679709813492295449586031792657667830635164741 / 24052915052186143128797603913750000000000000000000000000000000000\) | 是 |
| \(1/10000\) | \(23821519649033665711588756952800388139981734716562153495415164741 / 24052915052186143128797603913750000000000000000000000000000000000\) | 是 |

完整的 \(\kappa_{\max}\)、尾界和严格余量由脚本输出为准确分数；表中不以小数承担
严格比较。命令为：

```sh
TMPDIR=/Volumes/data/CodexTemp/T2 TMP=/Volumes/data/CodexTemp/T2 \
  PYTHONDONTWRITEBYTECODE=1 python3 -I -B \
  planar-centered-square-weak11/computation/check_supercell_mass_domination_tail.py
```

该命令同时运行一个非单位阈值回归：\(m=P_x=P_y=R=1\)、
\(\lambda=3\)、\(\delta=1/5\)、\(\kappa=7/8\)，其中
\(G_\delta=289/100\)、\(\kappa G_\delta<3\)。它防止把 \(\lambda\) 或
\(\kappa\) 从通用公式中遗漏。

## 5. 对后续计算的接口和边界

若候选通过 (9)，完整坐标盒枚举可在目标超胞元
\([0,s_xP_x)\times[0,s_yP_y)\) 和半径 cap \(R\) 上进行，随后用目标**实际**质量
\(m_\nu\) 形成比值 \(\lambda A_\nu^-/m_\nu\)；有限胞元截断和平滑也以该质量为
分母。本稿中的 \(\kappa\) 只是一条充足尾界，不能用 \(\kappa s_xs_ym\) 代替
\(m_\nu\)，更不能把 (9) 当成新的 \(C_2\) 界。

`vertical_compute` 可据此实现独立的 `periodic_halo_dominated.py`：逐标签检查基线
质量正性、目标质量非负、\(w'\le\kappa w\)、一致的 lift 后 \(\ell_\infty\) 位移
不超过 \(\delta\)，并报告目标实际质量、(9) 的准确余量和完整 halo 面积。旧的
输运证明、旧周期面积核及其冻结证书均不应就此改写。

