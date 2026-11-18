# 任意频率商的布尔切片、占用对偶与影响力

设 $F:X\to\{\pm1\}$，频谱 $S$，$s=|S|$，$d=\dim\operatorname{span} S$。任选物理子空间 $L\le X$ 与补坐标，将 $F$ 写作 $F(y,x)$，$x\in L$。对 $x$ 的部分变换系数为 $a_\beta(y)$，$B=\{\beta:a_\beta\not\equiv0\}=S|_L$，$p_\beta=\Pr(a_\beta(y)\ne0)$。真谱按 $\beta$ 分带，每带应用不确定性给
$$ s\ge\sum_{\beta\in B}\frac1{p_\beta}. $$
投影秩 $r=\dim\operatorname{span} B$ 满足 $d\le\operatorname{codim} L+r$。下述界不要求决策树、仿射分割或独立翻号。

## 1. 切片支撑集合的最优占用下界

令 $\mathcal F(B)$ 是所有包含于 $B$ 的真实布尔谱支撑族，或其任一指定外包；每个单点本身合法。定义
$$
\kappa(B)=\min_{p\in\operatorname{conv}\{1_T:T\in\mathcal F(B)\}}
\sum_{\beta\in B}p_\beta^{-1},
$$
零分量对应正无穷。真实的 $p$ 属于该凸包，故 $s\ge\kappa(B)$。其精确对偶为
$$
\sqrt{\kappa(B)}=
\max_{w\ge0:\ \sum_{\beta\in T}w_\beta\le1\ (T\in\mathcal F(B))}
\ \sum_{\beta\in B}\sqrt{w_\beta}.
$$
一向由 $\sum w_\beta p_\beta\le1$ 和 Cauchy 得到。反向取使目标最小的 $p$；因单点可混合，最小值有限且各分量正。沿 $1_T-p$ 的方向导数非负，说明 $\sum_{\beta\in T}1/p_\beta^2\le\sum_\beta1/p_\beta$。以 $w_\beta=(1/p_\beta^2)/(\sum_\gamma1/p_\gamma)$ 即达界。

若 $B$ 分为若干组且每个允许支撑只在一组内，取各组允许的子族，则 $\sqrt{\kappa(B)}=\sum_g\sqrt{\kappa(B_g)}$，由对偶约束独立直接得到。于是若 $r>C\sqrt s$，至少有一个组满足 $r_g>C\sqrt{\kappa(B_g)}$。这是凸占用模型的必要核；选中组的活动物理区域一般未必是仿射域，不能把它直接当成 $F$ 的真实布尔限制。

系数平方恒等式迫使任一切片支撑 $T$ 中每对不同点的差有第二个无序配对。若图 $\Gamma_B$ 将同属 $B$ 内某个四点平面的点相连，$T$ 必为其团。故可取图的连通分量作上述分组；给任何投影谱基 $I$，又可令 $w=1_I/\omega(\Gamma_B[I])$，得到 $r\le\sqrt{\omega(\Gamma_B[I])s}$。对只关心 $r$ 的问题，若另付 $h$ 个频率方向可令该团数至多 $k$，则总秩不超 $h+\sqrt{ks}$（这里 $h$ 包括所商去的全部方向）。

## 2. 基坐标的加权影响力

选 $B$ 的线性基 $\beta_1,\ldots,\beta_r$，把每个标签写在这组基下，坐标支撑记为 $J(\beta)$。切片只依赖这 $r$ 个二元变量。令
$$ m_i=|\{\beta:i\in J(\beta)\}|,\qquad
 d_i=\max_{i\in J(\beta)}|J(\beta)|,\qquad
 \gamma_i=\max\{1/m_i,2^{1-d_i}\}. $$
任何 $z_i\ge0$ 若满足 $\sum_{i\in J(\beta)}z_i\le1$ 对所有 $\beta$，则
$$
\sum_{i=1}^r\sqrt{z_i\gamma_i}\ \le\ \sqrt{s}.
$$
固定切片，相应 $\beta_i$ 系数活跃时，其第 $i$ 坐标奇部去除 $i$ 字符得到非零的 $\{0,\pm1\}$ 导数，剩余谱点至多 $m_i$，Fourier 次数至多 $d_i-1$。不确定性给该导数的非零概率至少 $1/m_i$。任意非零多线性多项式的度数 $k$ 给非零概率至少 $2^{-k}$：取最高次单项式，固定该单项式以外的全部变量，最高次系数固定非零，故在这 $k$ 个变量的至少一个取值处非零。这给另一项下界。非零概率等于影响力 $\operatorname{Inf}_i$。而
$$
\sum_i z_i\operatorname{Inf}_i=
\sum_\beta |a_\beta(y)|^2\sum_{i\in J(\beta)}z_i\le1 .
$$
平均得 $\sum_i z_i\gamma_i p_{\beta_i}\le1$，配合逐带不确定性和 Cauchy 即证。

一般可取 $z_i=1/d_i$。若所有投影标签次数至多 $k\ge1$，便有 $r\le\sqrt{k2^{k-1}s}$；更一般满足 $d_i\le k$ 的基坐标数至多这个数。该预算是给定商和基的充分证书，不保证任意谱都存在低次数商。

## 3. 四元符号约束能排除团图无法排除的切片

若某差在 $B$ 中恰有两组无序配对，则切片一旦含其中一对，便须包含全部四点，两个系数乘积之和为零。因系数为实数，四点的符号负号数为奇数。因而若某批被强制的全活四元组中，奇数个四元组指标异或为零，任何布尔切片都不可能同时触发它们。这给 $\mathcal F(B)$ 的有效外包收紧。

用独立方向 $e_{g,i},\lambda_g$ 构造标签包
$$
B_g=\{e_{g,i},e_{g,i}+\lambda_g:1\le i\le m_g\},\quad B=\bigcup_gB_g,\qquad m_g\ge2.
$$
任意跨包标签对在 $B$ 中的差只有一个配对，故切片只在一个包。包内若某切片同时使用至少三个不同基指标 $i,j,k$（每个指标先任选一个实际活跃标签），两指标所选点之差只有该对及另一平行对两种表示，所以两对的完整四边形都强制活跃。三个四边形的指标异或为零，三个符号奇校验矛盾。也可固定 $\lambda_g$ 的物理值：包函数在两个状态下都是独立 $e_{g,i}$ 字符的线性组合，布尔性迫使每个状态仅一个字符。故每一切片的支持至多四点，且要么单点，要么同包两指标的完整矩形。每个这样的矩形确实由 $\lambda_g$ 选择两个 $e$ 字符实现。

团图在单包内可有无界大团，而真实切片只能取上述四点。单包的凸占用值恰为 $m_g^2$：均匀混合所有指标对，使每个标签的活动率为 $2/m_g$；反向用每切片至多四点和 Cauchy。由第1节，$\sqrt s\ge\sum_g m_g$。投影秩为 $r=\sum_g(m_g+1)$，故若 $m_g\ge m\ge2$，则 $r\le(1+1/m)\sqrt s$。有全部 $m_g$ 偶且 $\sum m_g/2=2^h$ 的实际取等稀疏度例：各包指标配成不交对，用 $h$ 个外地址在所有对中选一对，在对内以相应 $\lambda_g$ 选择 $e_i$ 或 $e_j$ 字符。不同对的四个数据频点不交，每个点有 $2^h$ 个外地址频率，所以 $s=(\sum_g m_g)^2$。特别所有 $m_g=2$ 时，真实 $r=3\sqrt s/2$。以上不含商掉方向的维数，对总谱秩仍应加上该项。
