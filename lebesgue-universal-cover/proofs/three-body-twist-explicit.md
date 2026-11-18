# 三体零转角稳定性的显式角域

沿用零角切片的 $M_0$ 和三重接触局部证明的记号。五个参数按 $(a,y,\phi,b,z)$ 排列，角度取弧度。下面将局部证明量化到两个闭五维盒，并以外域内扇证书补齐全部平移。

## 两套解析接触方程

令 $p=(R_3+a,y)$、$U=a-R_3/2$、$u=(U,r+y)$、$\bar u=(U,-r+y)$，其中 $r=1/2$。设 $\rho=|p|$、$K=\sqrt{\rho^2-r^2}$、$e=p/\rho$、$f=Je$、$f_5=R_\phi(0,1)$，并记 $d=f_5\cdot f$, $o=f_5\cdot e$。两斜接触的法向和圆心可直接写为
\[
 n_\pm=(re\pm Kf)/\rho,\quad w_0=(-\rho+2gKd)e+(gr o/K)f,\quad w_\pm=w_0\mp gf_5.
\]
将表达式代入即可得到 $p\cdot n_\pm=r$ 及 $w_\pm\cdot n_\pm=-r$。因此五边形平移是 $t_5=w_0+R_\phi(C,0)$，其上下右顶点分别为 $v=w_0+R_\phi(C+k,r)$ 与 $\bar v=w_0+R_\phi(C+k,-r)$，其中 $g,C,k$ 的定义见零角局部证明。

第三条公共支持法向 $N$ 接近正竖直。参数化 $T=\arg N=\pi/2+2\arctan t$，其中
\[
 N=\frac{(-2t,1-t^2)}{1+t^2},\qquad t=\frac{y}{U-\sqrt{U^2+(y+m)y}}.
\]
这里 $m=1$ 对应 $u\cdot N=r$，而 $m=-1$ 对应 $\bar u\cdot N=-r$，即 $B_3$ 使用以 $\bar u$ 为心的单位圆弧。两式均由对应支持关系化作 $(y+m)t^2+2Ut-y=0$ 并选择近零根得到。最终用方程 $E_m(a,y,\phi)=v\cdot N-r=0$ 决定 $y$。两个分支在 $y=0$ 时方程值相同；严格正的 $E_{m,y}$ 说明二者的根总是同号。根为负时选 $m=1$，根为正时选 $m=-1$。在此选择下，分母严格负，$t$ 的符号确定 $T$ 位于 $\pi/2$ 的哪一侧。证书在两个解析分支上都证明 $|\partial_\phi T|<1$，由 $T(a,0)=\pi/2$ 得 $T<\pi/2+\phi$ 对所有正角成立，因而 $B_5$ 的上侧材料确实是顶点 $v$。

对每个有理 $(a,\phi)$ 子盒，把 $y$ 初步限在 $[-.001,.001]$，由 Arb 认证 $E_y>0$。取子盒中心求三步 Newton 近似 $\widehat y$，用 $(a,\phi)$ 两个一阶导数的区间笼住 $E(a,\widehat y,\phi)$。若其绝对值上界小于 $E_y$ 的下界乘 $\widehat y$ 到 $y$ 区间端点的距离下界，端点符号必异，根存在且唯一。再用均值定理收紧到 $y\in\widehat y-E(a,\widehat y,\phi)/E_y$。小盒之间的根由唯一性拼合。所有导数与面积系数用 160 位 Arb 级数核验，Newton 近似本身不承担证明。

## 锚点面积的显式二阶界

$m=1$ 的材料闭链使用零角局部证明的弧表。$m=-1$ 时作两处替换：在上侧 $T$ 到 $\pi/2$ 增加以 $\bar u$ 为心的单位弧，随后用原来的 $u$ 顶点；在下侧，$\bar u$ 顶点段直接延至 $\pi+T$，删去以 $u$ 为心的短弧。其余材料及连接同前文。程序核验各角的原有严格排序和对五边形扭角的适配。选择物理分支时新增短弧的方向正确。每个分支可以单独解析延拓到全角区间，下面的二阶估计在两个延拓上分别成立，因此选择分支时无需拼接导数。

为减小区间相关性误差，将 Green 面积简化如下。对支撑行 $Q=C_i+\rho_i n$，弧项为 $\rho_i\det(C_i,n_2-n_1)+\rho_i^2(\theta_2-\theta_1)$。所有异体连接都在共同 $D$ 支撑线 $n\cdot Q=r$ 上，同行材料衔接处连接项为零。故连接项可统一写为 $r(C_{i+1}-C_i)\cdot Jn$，总和把弧项的 $\rho_i$ 改为 $\rho_i+r$，不改变角长度项。再用 $\arg n_- -\arg n_+=-2\arctan(K/r)$、另一根方向 $2\arctan(w_{\pm,y}/w_{\pm,x})-\arg n_\pm$，得到倍面积
\[
\begin{split}
2H_m(a,\phi)={}&14\pi/15-2r^2\arctan(K/r)
 +2(1-r^2)\left[\arctan\frac{w_{-,y}}{w_{-,x}}-\arctan\frac{w_{+,y}}{w_{+,x}}\right]+\phi+I_m\\
 &+K+4g\rho d+4g^2 Kd^2-2g^2r o^2/K+2g(C+k)-\cos(\pi/6)-\cos(3\pi/10)\\
 &+2\det(v-U_m,N),
\end{split}
\]
其中 $(I_1,U_1)=(0,u)$、$(I_{-1},U_{-1})=(\pi-2T,\bar u)$。上式的圆弧角长度和来自三个盘弧合并、两个斜分别接点及上／下追加短弧。行列式合并中，$B_3$ 斜项为 $2\det(p,n_+-n_-)=4K$，$B_5$ 斜项为 $3[\det(w_+,n_+)-\det(w_-,n_-)]/2=-3K+6g\rho d$。其余非常数项来自顶部接点与 $2g(t_5\cdot R_\phi(1,0)+k)$；后一个表达式展开后得到第二行。常数由 $B_3$ 的 $150^\circ,210^\circ$ 弧端和 $B_5$ 的 $54^\circ,270^\circ,306^\circ$ 弧端给出。证书在每个子盒上按 $E_y y'+E_\phi=0$ 及二次项等于零恢复隐函数的二阶级数，再代入此面积式。

两个分支在零角均与实际凸包一致，零角所用的支撑顺序符号也在整个 $a$ 区间重验。因此 $H_m(a,0)\ge M_0$。两种分支的零角一阶项都为零，因为各自增加的短弧在零角退化为点，支撑点的切向滑动无一阶面积贡献，刚体速度在反射下为奇。对整个区间积分二阶下界即可认证 $H_m(a,\phi)\ge M_0+\kappa\phi^2$。宽盒可取 $\kappa=.11$。

## 平移盒内的定量下界

如局部证明，在六条三重线上各允许独立插入中间材料点。正面积角速度来自圆弧支持为正，所有异体连接和插点之后的连接端点的严格正行列式均在包含锚点与所有目标平移偏差的区间上核验。闭链绕数从零角凸包或同一闭链删去退化弧的极限保持为一，故每条链给出有效内扇面积。材料在两套物理分支的转换点衔接一致，因为增加的弧在该处变成点。三个代表连接上的中间材料次序由对应两段行列式正性保证，且三材料按接触方程共线。

记横截面偏差 $\xi=(\Delta y,\Delta b,\Delta z)$，不插点闭链的增量是线性形式 $\ell(\xi)$。其中潜在二次项为顶部和底部相反连接的两项 $\pm\Delta y\Delta b/2$，恰相消。三条代表线选择 $n_+,n_-,N$。给相反支持线一起比较，任选是否插入可获得增量 $|L_i(\xi)+Q_i(\xi)|$，因为相反线三个对应材料只差同一单位法向，且顺序反转，插点差对同一平移严格变号。$L_i$ 是插点面积增量的线性项，二次项依序为 $Q_i=(-1,1,-1)_i\,\Delta y\Delta b/2$。因此
\[
 A\ge H_m+\ell(\xi)+\sum_{i=1}^3|L_i(\xi)+Q_i(\xi)|.
\]

把 $L_i$ 的系数组成矩阵 $B$，用有理区间 Cramer 公式写 $\ell=\lambda^T B$。以 $w$ 包围 $\|\lambda\|_\infty$，以 $c_0$ 下界 $1/\max_{i,j}|(B^{-1})_{ij}|$，则 $\sum_i|L_i(\xi)|\ge c_0\|\xi\|_\infty$。若 $d_0$ 包围 $\min\{\sup|\Delta y|,\sup|\Delta b|\}$，则
\[
 A\ge H_m+\left[(1-w)(c_0-3d_0/2)-d_0|\lambda_2-\lambda_1-\lambda_3|/2\right]\|\xi\|_\infty.
\]
方括号已逐盒认证严格为正。这里实际考察的两个平移盒分别为：

| 模式 | $a$ | $y$ | $b$ | $z$ | $\phi$ |
| --- | --- | --- | --- | --- | --- |
| wide | $[-.017,-.0076]$ | $[-.0045,.0045]$ | $[.0185,.0274]$ | $[-.0045,.0045]$ | $[0,.0025]$ |
| near | $[-.0145,-.0094]$ | $[-.00076,.00076]$ | $[.0219,.02346]$ | $[-.00154,.000385]$ | $[0,.0045]$ |

## 角带的其余平移

第35轮对所有面积至多 $.838$ 的放置给出平移代表根盒。用第45轮冻结的完整五维 Tomar 种子树，并重新核验成员、绕数及面积，以 $0.834780945776$ 为严格下端。程序先将种子叶交角带，再去掉位于参数盒表中已证部分的内接 dyadic 盒。去除由四个平移坐标依次作切割，所得补集盒和内盒是完整角带的有限闭覆盖，另用精确有理体积求和验收。补集每个正体积盒使用 Taylor 检查或多扇有理权重检查，没有采纳原种子阈值所得的面积结论。

首个角带是 $[0,163/65536]$，补集 1666 个输入格、1682 次核验、0 残盒。第二个角带是 $[163/65536,147/32768]$，深度10首轮补集核验 4092 次，剩余59盒再经397次核验全部关闭。宽盒和窄盒分别包含两个角带所去除的 dyadic 盒。由此，对 $0\le\phi\le 147/32768$ 的所有平移，都有面积至少 $M_0$，正角严格成立。还能给出统一的二次增量 $10^{-8}\phi^2$，因为两个内盒的增量系数都大于 $0.11$，且 $10^{-8}(147/32768)^2<0.834780945776-0.8347809457752037$。对相反角度应用共同反射，得到显式双侧区间。第45轮最难父盒角域 $[0,2635359/1073741824]$ 完全落在此区间内；此次结论关于这段角域的全部平移。

复现需含 `python-flint` 的 Python。固定 SOURCE 为第39轮冻结 Tomar 检出，以下命令均从仓库根目录运行。

```sh
python -B lebesgue-universal-cover/computation/certify_three_body_twist_box.py wide --a-count 940 --phi-count 34
python -B lebesgue-universal-cover/computation/certify_three_body_twist_box.py near --a-count 510 --phi-count 64
python -B lebesgue-universal-cover/computation/certify_three_body_twist_exterior.py SOURCE wide OUTPUT_WIDE
python -B lebesgue-universal-cover/computation/certify_three_body_twist_exterior.py SOURCE near OUTPUT_NEAR
python -B lebesgue-universal-cover/computation/certify_three_body_twist_residual.py SOURCE OUTPUT_NEAR OUTPUT_CLOSURE
```

实际环境为 Python 3.13.18、python-flint 0.9.0、Arb 160 位。两个解析盒分别核验 63920、65280 个双分支单元。对应输出为 `three-body-twist-box-wide.json`、`three-body-twist-box-near.json`、`three-body-twist-exterior-wide.json`、`three-body-twist-exterior-near.json` 和 `three-body-twist-exterior-near-closure.json`。
