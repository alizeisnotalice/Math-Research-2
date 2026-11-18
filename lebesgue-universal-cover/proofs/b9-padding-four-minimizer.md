# 九边形相对四体极小凸包所需的外扩半径

设 $K_*$ 为圆盘与定宽三、五、七边形的四体极小凸包，取驻点证书中的规范代表。记 $\mathbb B$ 为单位半径的圆盘，定义
$$
 \rho_9(K)=\inf_{\theta\in\mathbb R,\ t\in\mathbb R^2}
 \inf\{\varepsilon\ge0:t+R_\theta B_9\subseteq K+\varepsilon\mathbb B\}.
$$
有严格夹逼
$$
 0.00394138015<\rho_9(K_*)<0.00394138020.
$$
当四体放置的凸包 $K$ 面积不超过 $0.836490455462$ 时，还可统一得到 $\rho_9(K)>0.00392411$。

## 全角下界

对三个有理单位方向 $u_i$，取满足 $\lambda_i>0$、$\sum_i\lambda_i=1$、$\sum_i\lambda_iu_i=0$ 的系数。系数由三个有向行列式归一化后得到，平移抵消用有理算术直接核对。若 $z_i\in R_{\theta_0}B_9$，则对任何平移 $t$ 和 $\theta=\theta_0+s$，有
$$
 \max_{|u|=1}\bigl(h_{t+R_\theta B_9}(u)-h_{K_*}(u)\bigr)
 \ \ge\
 \sum_i\lambda_i\bigl(\langle R_sz_i,u_i\rangle-h_{K_*}(u_i)\bigr).
$$
因此只需在每个角区间给出三个方向和各自的一点材料，而右式无需平移变量。

以 $2\pi/9$ 为旋转周期，按其有理倍数给出 41 个闭子区间。证书核对端点从周期的零端到单位端完全相接。各子区间中点为 $\theta_0$，半宽为 $w$。程序在中点选取每个方向上的顶点或圆弧材料；圆弧法向分区若无法定枝，就使用合法顶点。三个旋转贡献合并为
$$
 A\cos s+B\sin s,\qquad
 A=\sum_i\lambda_i\langle z_i,u_i\rangle,\quad
 B=\sum_i\lambda_i\det(z_i,u_i).
$$
所有子区间都有 $A>0$ 及 $0<w<\pi/2$，故在整个子区间可以用 $A\cos w-|B|\sin w$ 下界。这一估计保留不同方向一阶项的抵消。

宿主支持函数逐体取最大值。对奇数 $n$ 的正定宽 $n$ 边形，顶点 $v_k$ 的法向区间以其正多边形方位居中，圆弧中心 $v_k$ 对应的法向区间在相反方位居中，两类区间的半宽均为 $\pi/(2n)$。相应支持函数分别是 $\langle v_k,u\rangle$ 和 $1+\langle v_k,u\rangle$，其中顶点坐标已含该体的平移。程序枚举区间内所有可能的分支并求上界。所有圆盘方向的支持值均为 $1/2$。该全角核验使用 49 个有理单位方向，最小下界大于 $0.00394138015485$。

## 可行上界

取精确十进制放置参数
$$
 (\theta,t_x,t_y)=(0.35681038,\ 0.003790506253,\ -0.001905119134).
$$
用驻点参数区间核验宿主边界拓扑后，把其支撑链的切换法向与该放置九边形的切换法向求并，得到覆盖整个圆周的 37 个原子区间。每段的支持差具有形式 $a+\langle v,u\rangle$。该段最大值可在端点取得，或在 $v/|v|$ 落入端点所界的法向扇时取得。两种端点值逐一核验，行列式符号无法排除内部极大方向时就计入 $a+|v|$。核验的最大上界小于 $0.003941380194195$，所以该放置在规定上端的平行外扩内。

## 低面积四体凸包的稳定性

全域四体窄窗证书把面积不超过 $0.836490455462$ 的每个放置的一个规范代表送入定角平移盆地。对该代表，盆地的面积余量给出锚定法向的欧氏偏移 $<3.968\cdot10^{-6}$，及六个固定角平移偏差之和 $<4.373\cdot10^{-9}$。

锚定角参数到相应刚体的 Hausdorff Lipschitz 常数可由 $\|D t_n\|_{\mathrm F}+R_n\|\nabla\phi_n\|_2$ 上界，其中 $R_n=1/(2\cos(\pi/(2n)))$ 为含该中心化刚体的圆盘半径。圆弧法向区间内的支持值至多 $1/2$，顶点法向区间内的支持值至多 $R_n$，故该半径界成立。对连结两锚定法向的整条线段，用次数 8 的相关 Taylor 包络核验三个刚体的最大 Lipschitz 常数小于 $4.351$。取凸包后，Hausdorff 距离不超过各刚体的最大 Hausdorff 距离。所以
$$
 d_H(K,K_*)<4.351\,(3.968\cdot10^{-6})+4.373\cdot10^{-9}
 =0.000017269141
$$
对规范代表成立。更换宿主时，$\rho_9$ 的改变量至多为宿主的 Hausdorff 距离，这由单位方向支撑差逐项比较后对放置取下确界得到。将上述距离界从严格下端扣除，即得稳定版。九边形的镜面对称保证规范代表中的转动和平移仍覆盖原来允许的全部副本。

## 复算

在仓库根目录执行

    python -B lebesgue-universal-cover/computation/certify_b9_padding_four_minimizer.py \
      lebesgue-universal-cover/computation/b9-padding-four-minimizer.json

重放使用 Python 3.13.18、python-flint 0.9.0、Arb 160 位。冻结输入为 [b9-padding-four-minimizer.json](../computation/b9-padding-four-minimizer.json)，冻结输出为 [b9-padding-four-minimizer-audit.json](../computation/b9-padding-four-minimizer-audit.json)。数值线性规划仅用于寻找三个方向，严格重放只使用冻结后的有理方向和 Arb 区间。
