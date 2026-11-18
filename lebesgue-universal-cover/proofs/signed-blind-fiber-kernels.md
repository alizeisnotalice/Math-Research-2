# 盲叶 1202/32259 沿五边形横移分四类后的共有内核

设 \(K(q)\) 为圆盘和四个正 Reuleaux 多边形在原叶五体位姿 \(q\) 的凸包，\(T=0.83743\)。若 \(K(q)\) 已有面积至少 \(T\)，加入任意第六个输入体也不能降低面积。以下针对第 1202 分片的前序节点 32259，以有符号面积下端处理 \(A(K(q))<T\) 的必要条件，再检验按一个平移坐标分成四类是否足以使用少数共有内核。

## 五体筛选与共有内核

以冻结原树恢复闭盒。取其名义材料链上八个平移项的一阶系数绝对值乘相应坐标宽度，最大项为 \(B_5\) 的横移 \(t_{5,x}\)。按该坐标均分闭区间，相邻类别共享端点。其五个整数端点为
\[
390479953920,\quad 404425666560,\quad 418371379200,\quad
432317091840,\quad 446262804480
\]
除以 \(2^{44}\)。沿用五体面积的平移端点判定，将每类余盒分支至深度 11。比较整数目标 \(\lfloor 2^{57}T\rfloor+1\)，通过的盒全盒满足 \(A(K(q))>T\)。任何五体面积小于 \(T\) 的位姿一定在其类别的至少一个保留闭子盒内。

| 类别 | 调用次数 | 关闭盒 | 保留盒 | 仅由有符号改进关闭 |
|---:|---:|---:|---:|---:|
| 0 | 167 | 66 | 18 | 10 |
| 1 | 369 | 122 | 63 | 13 |
| 2 | 721 | 210 | 151 | 33 |
| 3 | 1069 | 293 | 242 | 35 |

对每个保留盒 \(B\)，在其平移中点放置冻结材料表中的每个体。表行的名义转角取最靠近盒角中点的可用表行。设 \(\epsilon\) 为真实平移距该中点的欧氏距离上界，加上真实转角与表行之差乘该体支撑函数的角 Lipschitz 上界。设 \(rB_2\) 包含于该名义体。若名义体的支撑为 \(h\)，扣去平移中点后，真实体的支撑至少为 \(h-\epsilon\)，而
\[
(1-\epsilon/r)h\le h-\epsilon
\]
在单位方向成立。将名义材料点以 \(1-\epsilon/r\) 关于原点缩放，再平移到盒中点，所得点因而落在真实体中。四个 Reuleaux 体用 \(r=0.42,0.47,0.48,0.49\)，三个可旋转体分别用角 Lipschitz 上界 \(0.163,0.115,0.089\)。这些界由 \(1-\frac1{2\cos(\pi/(2n))}>r\) 与 \(\frac12\tan(\pi/(2n))<L\) 认证。以固定圆盘的材料点和各缩入体材料点的有理凸包构成 \(R_B\)。因此 \(R_B\subset K(q)\) 对全部 \(q\in B\) 成立。

每类共有内核 \(P_i\) 取该类别全部 \(R_B\) 的交。由于各 \(R_B\) 均包含同一个原点在内部的圆盘采样多边形，程序以有理凸包计算 \(\left(\operatorname{conv}\bigcup_B R_B^\circ\right)^\circ\)，得到 49、54、55、54 边。其有理面积的十二位小数下端依次为 \(0.829800199521,\ 0.829152166325,\ 0.827600259599,\ 0.827042992054\)。

## 平衡共享容量的上端

固定一个类别的内核 \(P=P_i\)。以 \(n_j\) 表示外边法向，其长度为对应边长。对平衡权重行 \(\lambda_{aj}\ge0\)，有 \(\sum_j\lambda_{aj}n_j=0\)，每条边允许总权重 \(\sum_a\lambda_{aj}\le1/2\)。在任意给定角度处，给每个测试体任选平移 \(s_a\)，并设
\[
\widetilde H=\operatorname{conv}\!\left(P\ \cup\ \bigcup_a (R_aX_a+s_a)\right).
\]
行的值不随 \(s_a\) 改变，因为 \(\sum_j\lambda_{aj}\langle s_a,n_j\rangle=0\)。逐个方向应用 \(\widetilde H\) 的支撑上端及 \(\widetilde H\supset P\)，给出
\[
|P|+\sum_{a,j}\lambda_{aj}\big(h_{R_aX_a}(n_j)-h_P(n_j)\big)
\le |P|+\tfrac12\sum_j\big(h_{\widetilde H}(n_j)-h_P(n_j)\big)
= V(P,\widetilde H).
\]
人工平移仅用于这个上端。五体转角由同类别的真实原盒点固定，第六体用冻结的有理 P30 点集和一个明确角度。P30 关于横轴反射不变，两手性共用同一允许角集合。对整类有效的角分段必须包含所选角度，因此任意加密分段和权重行仍受该上端约束。每条边以 192 位 Arb 认证支撑并向上取整，最终上端如下。

| 类别 | 见证的真实五体面积上端 | \(V(P_i,\widetilde H)\) 上端 | 低于 \(T\) 至少 |
|---:|---:|---:|---:|
| 0 | 0.837540228314 | 0.836617816077 | 0.000812183923 |
| 1 | 0.837387365851 | 0.836394447929 | 0.001035552071 |
| 2 | 0.837336262352 | 0.836119720687 | 0.001310279313 |
| 3 | 0.837325631020 | 0.836109464855 | 0.001320535145 |

类别 0 的见证在有符号筛后超集，其真实五体面积是否小于 \(T\) 未获认证。类别 1、2、3 的见证通过 \(k=1024\) 的真实五体上端验证，均严格在五体低面积子集。因此这个四类共有内核构造无法用平衡混合面积证书替代盲叶的第六体树，即使在后三类中只要求覆盖真实低面积子集。结论适用于这四个明确内核。五体筛选后允许另选任意共有内多边形的最优预算，以及 \(t_{5,x}\) 任意细分后保持同一策略的能力，仍未由此界定。

## 复算

三个入口只进行原树整数筛选、有理凸包、支撑区间和真实五体面积上端。选择见证时使用的数值优化不参与验证。原树采用冻结的 `cert11d_08366`，材料表采用 R87 的 `table32.npz`。三个输出与归档逐字节比对。工作目录 `WORK` 与编译模块缓存 `CACHE` 可位于 CodexTemp。

    mkdir -p WORK
    python -I -B computation/check_signed_fiber_classes.py COMPUTATION CERT TABLE WORK CACHE > WORK/scan.json
    cmp WORK/scan.json computation/signed-blind-fiber-classes.json
    python -I -B computation/construct_signed_fiber_kernels.py COMPUTATION CERT TABLE WORK/scan.json > WORK/cores.json
    cmp WORK/cores.json computation/signed-blind-fiber-kernels.json
    python -I -B computation/certify_signed_kernel_obstruction.py COMPUTATION WORK/scan.json WORK/cores.json computation/signed-blind-fiber-witnesses.json > WORK/caps.json
    cmp WORK/caps.json computation/signed-blind-fiber-obstruction.json

复算使用 Python 3.13.18、python-flint 0.9.0、NumPy 2.5.5、Apple Clang 21.0.0。原树 SHA256SUMS 的 SHA-256 为 `56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410`，R87 材料表为 `fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c`。
