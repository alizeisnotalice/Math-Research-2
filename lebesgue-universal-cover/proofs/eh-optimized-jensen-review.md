# EH 新弦的 Jensen 核验

取 \(\sigma=1.5494009336023418793^\circ\)。独立的 160 位 Arb 计算由支撑线和圆相交构造 \(U,T,J,C_3,B_3,S\)，用原第47轮面积上端作基线。原有五种放置论证沿用已建立的几何，引入新弦所需的方向和位置条件逐项用 Arb 核验。

## 两次 Jensen 与投影符号

设 \(K,Q_2K\) 各命中新耳一点，\(y,y'\) 是二点低于 \(U_1\) 的高度，\(r,r'\) 是各体在竖直边 CB 上的接触点高于 \(C_3\) 的高度。弦斜率用 \(k=(T_x-U_{1,x})/(U_{1,y}-T_y)>0\) 表示。设 \(p=u(7\pi/6)\)，\(\eta=1/2-p\cdot U_1\)，\(\gamma=1/2-\sqrt3 k/2\)。由 \(Q_2p=-p\) 和单位直径，\((y+y')/2\le\eta/\gamma\)。设 \(Y\ge\eta/\gamma\)，\(H=1/2-U_{1,x}\)，\(a=C_{3,y}-U_{1,y}\)，\(c=Hk>a>0\)。原证明的函数 \(f(y)=\sqrt{a^2+2cy+y^2}-a-y\) 满足

\[
f'(y)=\frac{y+c}{\sqrt{a^2+2cy+y^2}}-1>0,\qquad
f''(y)=\frac{a^2-c^2}{(a^2+2cy+y^2)^{3/2}}<0.
\]

由 \(r\le f(y)\) 和 \(r'\le f(y')\)，可取 \(R_0\ge f(Y)\) 作平均上界，\(R_{\rm ind}\ge f(2Y)\) 作逐点上界。选择 \(0<R_0\le R_{\rm ind}\)。

令 \(m=m_5\)，\(M=-m_y>0\)，\(b=C_3+r e_y\)，\(w=C_3-m\)。若 \(w_K\) 是 \(-m\) 支撑点，则 \(e=b-w_K-m\) 满足 \(-Mr\le e\cdot m\le 0\)。因为 \(MR_{\rm ind}<1\) 且 \(|b-w_K|\le1\)，有 \(|e\cdot J_0m|\le\sqrt{2Mr-M^2r^2}\)。对所选 \(v\)，保留 \(v\cdot m<0\) 的符号给出 \((v\cdot m)(e\cdot m)\le-Mr(v\cdot m)\)。以 \(\beta=v\cdot J_0m>0\)，\(\lambda=v_y+M(v\cdot m)>0\)，得到

\[
v\cdot(w_K-w)\ge-L(r),\quad L(r)=\beta\sqrt{2Mr-M^2r^2}-\lambda r.
\]

\(L\) 在 \(0\) 连续，对正 \(r\le R_{\rm ind}\) 有

\[
L''(r)=-\frac{\beta M^2}{(2Mr-M^2r^2)^{3/2}}<0,\qquad
L'(r)\ge\frac{\beta M(1-MR_{\rm ind})}{\sqrt{2MR_{\rm ind}}}-\lambda>0.
\]

故两点的平均损失不超过 \(L(R_0)\)。取 \(K\) 的 \(n_5\) 支撑点 \(N\)。\(Q_2\) 关于 \(n_5\) 轴反射，故 \((N+Q_2N)/2=q_0=n_5/2\)。将 \(Q_2K\) 中的第二个支持点及测量方向反射回 \(K\)，两个真实投影各至多一，平均至少

\[
v\cdot(w-q_0)-L(R_0)=v\cdot\xi,\quad
\xi=w-q_0+R_0e_y+MR_0m-\sqrt{2MR_0-M^2R_0^2}\,J_0m.
\]

选 \(v=\xi/|\xi|\) 并检验所需符号。在 \(\delta=0.0000154\) 时 Arb 得 \(v\cdot m<-0.92609\)、\(\beta>0.37729\)、\(\lambda>0.33172\)、\(L'>11.08324\)、\(|\xi|-1>4.78971354\times10^{-7}\)。故不可能重复命中。附加帽面积超过 \(3.4455208147077624\times10^{-9}\)，得到 \(\Lambda<0.844093590940968945\)。

## 保留未展开的平方距离

弦条件给 \(z_x-U_{1,x}\le ky\)，证书检验 \(H>2kY\)。于是

\[
1\ge|b-z|^2\ge(H-ky)^2+(a+r+y)^2,
\qquad r\le f_*(y)=\sqrt{a^2+2cy-k^2y^2}-a-y.
\]

被开方式在 \([0,2Y]\) 为正，且严格凹，故 \(f_*\) 严格凹。检验 \(f_*'(2Y)>0\) 即知整个区间递增。用 \(f_*\) 替换 \(f\)，两次 Jensen 及其余论证沿用。在 \(\delta=0.0000154059\) 时，\(|\xi|-1>7.25877398\times10^{-9}\)，附加帽面积超过 \(3.4468410320698719\times10^{-9}\)，得到 \(\Lambda<0.844093590939648727\)。\(\delta=0.000015406\) 时该投影判据的余量为负。

复现使用 Python 3.13.18 与 python-flint 0.9.0。由 `computation/audit_optimized_eh_chord.py OUTPUT --old-area computation/047-center.json` 重算四组试点。实际输出在 `computation/067-optimized.json`。
