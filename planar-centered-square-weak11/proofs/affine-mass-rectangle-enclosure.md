# 仿射守恒质量矩形族的完整 halo 外包络

日期：2026-11-15（Asia/Shanghai）。本稿给出固定原子坐标与周期下，一个二维仿射、
总质量守恒的参数矩形族的全局上包接口。它认证的是该**受限质量族**上的周期 halo
比值上界；不构造把各盒逐点最大值同时达到的单一输入，也不声称二维一般弱 \((1,1)\)
常数的上界或精确值。

## 1. 对象、量词与完整 halo 表示

固定 \(P_x,P_y,\lambda,R>0\)，基胞 \(F=[0,P_x)\times[0,P_y)\)，以及固定原子坐标
与固定周期副本窗口。令 \(\mathscr B\) 是完整坐标盒定理在半径 \(R\) 下的**所有**
非空源盒；写
\[
 B=[x_B^-,x_B^+]\times[y_B^-,y_B^+],\qquad
 M_B(\theta)=c_B\mathbin\cdot w(\theta),\quad \theta=(a,b). \tag{1}
\]
这里 \(c_B\) 是该盒中每个标签的固定非负重数。参数域为闭矩形
\(D=[a_-,a_+]\times[b_-,b_+]\)，质量为仿射函数
\[
 w_i(a,b)=w_i^0+u_i a+v_i b>0\quad(\theta\in D),\qquad
 \sum_iw_i(a,b)=M_0>0. \tag{2}
\]
假设已对**同一** \(D\) 认证统一尾界：对每个 \(\theta\in D\)，半径 \(r\ge R\) 不存在
严格 \(\lambda\)-见证。该假设不能由参数网格抽样代替。

令
\[
 t_B(\theta)=\min\!\left\{R,\sqrt{M_B(\theta)/(4\lambda)}\right\},\qquad
 E_B(t)=F\cap(x_B^+-t,x_B^-+t)\times(y_B^+-t,y_B^-+t). \tag{3}
\]
端点次序错误时 \(E_B(t)=\varnothing\)。完整坐标盒定理逐参数给出
\[
 H_\theta\cap F=\bigcup_{B\in\mathscr B}E_B(t_B(\theta)),\qquad
 J(\theta)=\frac{\lambda |H_\theta\cap F|}{M_0}. \tag{4}
\]
每个 \(M_B\) 是一个**共同**参数 \(\theta\) 的线性式；以后对不同 \(B\) 取得的最大值
绝不被解释为某个真实 \(w(\theta)\)。

## 2. 单参数 cell 的严格外包络

令 \(C\subset D\) 是闭参数小矩形，\(V(C)\) 是其四个角。仿射性给出
\[
 M_B^+(C)=\max_{\theta\in C}M_B(\theta)=\max_{v\in V(C)}M_B(v),\qquad
 t_B^+(C)=\min\!\left\{R,\sqrt{M_B^+(C)/(4\lambda)}\right\}. \tag{5}
\]
因此，对每个固定 \(\theta\in C\)，\(t_B(\theta)\le t_B^+(C)\)，而 \(E_B(t)\) 对 \(t\)
单调不减。定义
\[
 U_C=\left|\bigcup_{B\in\mathscr B}E_B(t_B^+(C))\right|. \tag{6}
\]
则严格地
\[
 \boxed{\quad J(\theta)\le \frac{\lambda U_C}{M_0}\quad(\theta\in C).\quad} \tag{7}
\]
证明只需把 (4) 中每一个接收矩形包含进 (6) 中同一 \(B\) 的矩形，再取并与面积。
分母是 (2) 的实际固定总质量 \(M_0\)，并非 \(\sum_i\max_Cw_i\)、\(\sum_BM_B^+\)
或任何外包络“质量”。特别地，(6) 是多个可能互不相同角点的矩形的并；它通常不是
任何 \(\theta\) 的 halo，也不需要是。

为避免数值根号的向外扩张把整 cell 上本来恒空的盒错误激活，先作如下**精确**剪除。
若
\[
 x_B^+-x_B^-\ge2t_B^+(C)\quad\hbox{或}\quad
 y_B^+-y_B^-\ge2t_B^+(C), \tag{8}
\]
则对所有 \(\theta\in C\)，\(E_B(t_B(\theta))=\varnothing\)。这种盒可在该 cell 删除；
记余下集合为 \(\mathscr B_C\)。等号也可删除，因为接收矩形使用严格开区间。这个剪除
只省略全 cell 的真实空矩形，故不改变 (7)，却不能用一个抽样点的空性替代。有理输入中
这一步也无需近似根号：cap 支直接比较跨度与 \(2R\)，未 cap 支将比较平方为
\(\lambda(x_B^+-x_B^-)^2\ge M_B^+(C)\) 或对应 y 式。

## 3. 有理 dyadic 向外根界与准确并面积

现进一步假定固定几何、\(\lambda,R\)、仿射系数及 cell 角均为有理数。对 \(k\ge0\)，令
\(D_k=2^k\)，并定义
\[
 q_{B,C,k}=D_k^{-1}\left(\left\lfloor D_k\sqrt{M_B^+(C)/(4\lambda)}\right\rfloor+1\right),
 \qquad\widehat t_{B,C,k}=\min\{R,q_{B,C,k}\}. \tag{9}
\]
于是
\[
 0\le\widehat t_{B,C,k}-t_B^+(C)\le D_k^{-1}. \tag{10}
\]
这里始终加一格，即使根恰为 dyadic 也向外扩张；故未 cap 时误差严格为正。若
\(M_B^+/(4\lambda)=A/B\)（既约正有理数），整数 \(\ell=\lfloor D_k\sqrt{A/B}\rfloor\)
由 \(\ell^2B\le AD_k^2<(\ell+1)^2B\) 精确核验；不可用向下舍入浮点根号。令
\[
 U_{C,k}=\left|\bigcup_{B\in\mathscr B_C}E_B(\widehat t_{B,C,k})\right|. \tag{11}
\]
所有端点现为有理数。将所有 x、y 端点与 \(0,P_x,0,P_y\) 排序，以相邻端点开小格
测试是否被至少一个闭包矩形覆盖，即可用有理算术精确求出 (11)；开闭端点的差为零面积。
由 (7)–(11)，
\[
 J(\theta)\le\lambda U_{C,k}/M_0\quad(\theta\in C). \tag{12}
\]

## 4. 网格收敛、全局上界和达到

设 \(\mathscr P_h\) 是将 \(D\) 分割为 \(\ell_\infty\)-直径至多 \(h\) 的有限闭矩形 cell
的任意网格。每个保留源盒都含至少一个标签；由 (2) 的严格正性和紧性，
\(m_B:=\min_D M_B>0\)。若
\(M_B(a,b)=\alpha_B+\beta_Ba+\gamma_Bb\)，则对 \(\theta,\theta'\) 同属一个 cell，
\[
 |t_B(\theta)-t_B(\theta')|
 \le \frac{(|\beta_B|+|\gamma_B|)h}{4\sqrt{\lambda m_B}}. \tag{13}
\]
这是未 cap 根式的均值定理与 \(t\mapsto\min(R,t)\) 的 1-Lipschitz 性；cap 切换不另生
误差。记右端对 \(B\) 的最大值为 \(\delta(h)\)，并令
\(\varepsilon_{h,k}=\delta(h)+2^{-k}\)。

若一个被 \(F\) 裁剪的接收矩形各方向半径增加至多 \(\varepsilon\)，其新增面积至多
\[
 2\varepsilon(P_x+P_y)+4\varepsilon^2. \tag{14}
\]
这是把宽、高分别至多增加 \(2\varepsilon\) 的粗上界。令 \(N=|\mathscr B|\)。固定
\(C\in\mathscr P_h\) 并任选 \(\theta_C\in C\)，逐矩形包含、再以并的次可加性得到
\[
 0\le U_{C,k}-|H_{\theta_C}\cap F|
 \le N\{2\varepsilon_{h,k}(P_x+P_y)+4\varepsilon_{h,k}^2\}. \tag{15}
\]
故由 \(\overline J_{h,k}=\max_{C\in\mathscr P_h}\lambda U_{C,k}/M_0\)，
\[
 \max_{\theta\in D}J(\theta)
 \le\overline J_{h,k}
 \le\max_{\theta\in D}J(\theta)+
 \frac{\lambda N}{M_0}\{2\varepsilon_{h,k}(P_x+P_y)+4\varepsilon_{h,k}^2\}. \tag{16}
\]
因 \(\varepsilon_{h,k}\to0\) 当 \(h\to0,k\to\infty\)，分支包络从上方一致收敛。
这里用任意 \(\theta_C\)，而非为每个盒拼造一个“最大质量输入”，正是 (16) 成立且不犯
共同输入错误的原因。

若实现固定使用 \(D_k=2^{40}\)，即使 \(h\to0\)，(16) 仍保留由 \(2^{-40}\) 引起的
显式误差地板；此时只能报告该精度预算内的认证上界，不能声称数值 gap 可任意小。只有
同时令 bits \(k\to\infty\)（或预先选择足以满足目标误差预算的 bits）才有上段的一致
收敛结论。

另一方面，(3) 的半径连续，有限矩形并面积随端点连续，故 \(J\) 在紧集 \(D\) 上连续。
因此 \(\max_DJ\) 确实达到；(16) 提供其可认证全局上界序列。它不提供解析闭式最大点，
也不把逐 cell 的上界误报为已达到的值。

## 5. 精确小例：不同盒的角点最大值并非真实输入

此例只核验外包络算术，不是一个项目候选的完整 halo 计算。取
\(D=[-1,1]^2\)、总质量 \(M_0=6\)、\(\lambda=1\)，两盒质量
\[
 M_1(a,b)=3+a+b,\qquad M_2(a,b)=3-a-b. \tag{17}
\]
它们均正且总和恒为 6；但 \(M_1^+=M_2^+=5\) 分别在 \((1,1)\) 与 \((-1,-1)\)
达到，不能同时由一个参数点达到。令 \(F=[0,4]^2\)，两退化源盒中心为 \((1,1)\)、
\((3,1)\)，并取 \(R=5\)。在 dyadic 精度 \(k=2\)，
\[
 \sqrt5/2<5/4,\qquad\widehat t_1=\widehat t_2=5/4.
\]
两个裁剪接收矩形的准确并面积为 \(U=9\)：它们各有面积 \(81/16\)，交面积为 \(9/8\)。
所以该单 cell 的严格认证上界为 \(\lambda U/M_0=3/2\)。
`computation/check_affine_mass_enclosure.py` 用 `Fraction` 核对四角最大值、dyadic 向外根界、
准确并面积和固定分母；它特意断言不存在同时令两盒质量等于 5 的参数点。

## 6. R20 接口与范围

R20 的候选输入为 \(a\in[-1,1]\)、\(b\in[0,10]\) 的 204 标签守质量仿射族。新 kernel
应逐 cell 输出：cell 四角、每个完整源盒的 \((\alpha_B,\beta_B,\gamma_B)\)、四角
\(M_B^+\)、dyadic \(\widehat t_{B,C,k}\)、准确 \(U_{C,k}\)、固定 \(M_0\) 与
\(\lambda U_{C,k}/M_0\)。在采用本稿前还须逐项核验：所有标签在整个参数矩形严格正、
标签质量和确为常数、\(\mathscr B\) 是完整盒集，以及 R8 的 \(\kappa\le1001/1000\)、固定
位移输运确实对全参数矩形给出同一 \(R=5\) 尾界。任一项仅在抽样点成立都不足以使用
(7) 或 (16)。
