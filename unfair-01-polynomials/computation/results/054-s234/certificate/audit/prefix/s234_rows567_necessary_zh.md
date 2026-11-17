# S234 的第 5--7 行：整数平移与 Newton 必要筛

本稿只给出七次实二次数域 primitive core 的必要条件和有限前缀筛，**不**声称 S234 或原均匀和猜想已排除。设

\[
 P=1+\ell x+a x^2+(\ell a+k)x^3+d x^4+e x^5+f x^6+x^7,
\]

其中 \(a,d\in\mathcal O_K\setminus\mathbb Q\)，\(\ell,k,e,f\in\mathbb Z\)，横线为另一实嵌入；此处沿用 S234 已建立的输入分支 \(\ell\ne0\)（\(\ell=0\) 不在本筛的输入分类中）。令
\(A=GP,B=Q\bar P\)，其中 \(G,Q\in\mathbb Z[x]\)、常数项为一，\(A,B\ge0\)，且 \(C=AB\) 二元。采用此前 S234 的低端归一化

\[
G_1=Q_1=-\ell,\quad
G_2=\ell^2+u_2,\quad Q_2=\ell^2+v_2,\quad
G_3=-\ell G_2-k,\quad Q_3=-\ell Q_2-k,
\]

及

\[
G_4=u_4+\ell^2G_2+2\ell k,\qquad
Q_4=v_4+\ell^2Q_2+2\ell k,
\]

其中所有 \(u_i,v_i\) 为整数。已知第 2、4 行给出
\(0<A_2,B_2<1\)、\(A_4,B_4\ge0\)，以及
\(A_2B_4+A_4B_2>0\)。

## 平移参数

定义整数 \(\eta,\rho,\gamma\) 使得

\[
\begin{aligned}
G_5&=\eta-e-kG_2-\ell G_4,\\
G_6&=\rho-f+2\ell e-\ell\eta+2\ell kG_2+\ell^2G_4+k^2.
\end{aligned}\tag{1}
\]

直接卷积给出

\[
\begin{aligned}
A_5&=\eta-ka-\ell d,\\
A_6&=a(u_4+\ell k)+d(u_2+\ell^2)+\rho.\tag{2}
\end{aligned}
\]

第 7 行可把全部有理项（包括自由整数 \(G_7\)）吸入某个整数 \(\gamma\)，从而

\[
A_7=a(\eta-e-kG_2)-d(\ell G_2+k)+\gamma.\tag{3}
\]

共轭侧有独立整数 \(\zeta,\sigma,\delta\)。这些替换只是整系数的一一重参数化；特别不能把 \(e\) 或 \(f\) 从母式本身删除。

第 5 行满足 \(C_5=A_5+B_5\in\{0,1\}\)。若其无理部分消失，闭单位区间的三个端点组合都必须保留；否则两项是严格互补。第 6 行为

\[
C_6=A_6+B_6+A_2B_4+A_4B_2=1,
\]

所以 \(A_6,B_6\in[0,1)\)，并且给定第 4 行后 \(\rho,\sigma\) 唯一。第 7 行为

\[
C_7=A_7+B_7+A_2B_5+A_5B_2\in\{0,1\};\tag{4}
\]

最后交叉项非负但不必严格正，故此行仍要枚举闭端点，而不能套用第 6 行的唯一性。

## Newton 界和完整整数盒

对倒数根的 Newton 递推独立展开为

\[
\begin{aligned}
s_5={}&5ka+5\ell d-\ell^5-5\ell^2k-5e,\\
s_6={}&-2a^3-6\ell ka+6ad-6\ell^2d+\ell^6+6\ell^3k+3k^2+6\ell e-6f,\\
s_7={}&-7ka^2+7\ell^2ka-7\ell ad+7ae-\ell^7-7\ell^4k+7\ell^3d-7\ell^2e-7\ell k^2+7\ell f+7kd-7.
\end{aligned}\tag{5}
\]

根环极点界给出

\[
|s_5|<\frac{3107}{32},\quad |s_6|<\frac{12355}{64},\quad |s_7|<\frac{49283}{128},
\]

且系数包络给 \(|e|<117/4\)、\(|f|<17/2\)、\(|\ell|<17/2\)，即

\[
e\in[-29,29]\cap\mathbb Z,\qquad f\in[-8,8]\cap\mathbb Z,\qquad \ell\in[-8,8]\cap\mathbb Z.\tag{6}
\]

因此在尚未枚举 \(e\) 前，\(s_5+5e\) 必须满足双嵌入严格三角界

\[
|s_5+5e|<\frac{3107}{32}+145=\frac{7747}{32}.\tag{7}
\]

固定 \(e\) 后，令 \(S_6=s_6+6f\) 和 \(S_7=s_7-7\ell f\)。必要的廉价筛为

\[
|S_6|<\frac{12355}{64}+48=\frac{15427}{64},\qquad
|S_7|<\frac{49283+7168|\ell|}{128}.\tag{8}
\]

随后对 (6) 的 **全部** 59 个 \(e\) 和 17 个 \(f\) 再直接检查 (5) 的三条严格双嵌入不等式。故 (7)--(8) 只用于加速，不替代完整 \(e,f\) 盒。

## 计算状态的正确计数

一个输出记录固定到第 7 行的

\[
(\ell,a,k,d,u_2,v_2,u_4,v_4,\rho,\sigma,\eta,\zeta,e,\gamma,\delta)
\]

并附带一个列表 `f_values`。该列表中的元素只是该同一前缀的可能六次系数；它们不是不同的第 7 行状态，更不是已经构成的不同母式。完整母式还须通过后续根环和任意余因子条件。

`verify_formulas_independent.py` 用 SymPy 1.14 实际重展开了 (2)、(3)、(5)，而 `s234_rows567_independent.py` 用自写的 \(\mathcal O_K\) 整数算术、无浮点符号判定和每域原子 JSON 输出执行上述必要筛。
