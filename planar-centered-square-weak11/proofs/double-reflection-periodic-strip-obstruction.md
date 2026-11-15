# R54：双反射周期条带的完整边界系数与有限 L¹ 传递

日期：2026-11-15。对象固定为二维 Lebesgue 接收、轴平行中心共同半径正方形、全部 r>0、严格阈值 λ=1。本稿冻结 R46，新增 R41 规范相位的左边界计算。这里的 J 始终是阈值 1 的面积/质量比，不是对同一输入再取 sup_λ 的弱范数。

## 1. 结论与范围

取冻结 R41 的 204 个胞元原子 (α_t,β_t,w_t)，原周期 Px、Py、质量 m，完整周期 halo 每胞元面积 A。对整数 N≥7 定义无限纵周期条带

    σ_N = Σ_{i=0}^{N−1} Σ_{j∈Z} Σ_t w_t δ_(α_t+iPx, β_t+jPy).

记完整条带 halo 真正左右端点为 ℓ、r_N=(N−1)Px+C，W_N=r_N−ℓ。关于 x=c 的反射为 T_c。双反射闭合定义为

    Π_N = Σ_{k∈Z} τ_(2kW_N) [σ_N+(T_(r_N))#σ_N].

它的横周期准确是 2W_N，纵周期 Py，每胞元质量准确是 2Nm。新的严格结论为

    J(Π_N;1) = A/m + K/(mN) < A/m − 12/(5mN),    N≥7.             (1)

其中 K=b_x+γ_L+γ_R，γ 是对应接缝半窗口的完整 new-minus-old 面积。准确有理账本给

    γ_L ∈ [1.0328056783975574, 1.032805679474352] （区间端点近似显示），
    K ∈ [−2.436023761291603, −2.4360237578629964] （区间端点近似显示），
    K < −12/5 （准确有理判断）。

所列小数仅用于显示；严格符号由结果 JSON 中的 Fraction 比较承担。每个固定 N 的规范有限截断并缩小平滑尺度的 L¹ 序列，其阈值 1 比值趋于 J(Π_N;1)，且最终小于 A/m−6/(5mN)。这不是所有有限截断、任意平滑尺度或所有阈值的排除。

本轮仅对 R41 做新左缝计算，没有把结论冒充 R35 新计算。N≤6、其他相位、其他周期模板及一般 C₂ 精确值仍未解决。

## 2. 不取模的左轴变换

令 X=max_t α_t=979821360288625294630001/150000000000000000000000。原输入 min α=0、X<Px。使用等距变换

    (x,y) ↦ (x′,y′)=(X−x,y).

所有新标签 α′=X−α 仍在 [0,Px)，且

    X−(α+iPx)=(X−α)+(−i)Px.

所以左半平面源 i≥0 准确变成新标签的右半平面源 i′≤0，没有取模、没有标签相位改动。原左半接缝窗口 [ℓ,ℓ+15]×[0,Py] 准确映为 [C′−15,C′]×[0,Py]，C′=X−ℓ。因此 frozen `reflection_edge_seam.py` 与 `reflection_edge_states.py` 的既有右半平面接口可直接复用。

变换后的输入 SHA-256 为 `deb65c6011d01f61253cb1db6b4b0e919f014ad32711677b9ee3fd80645694bf`，原输入 SHA 为 `ab714c7f457ccf69f68fd629c947021faab16631c8995caa16b67c1ca2752db3`。独立逐点核验保留全部 204 个质量、y 坐标与周期。不能把该结论替换为“先 x→−x，再逐标签取模并统一 cutoff”；后者一般错误，独审已有准确反例。

轴小核从完整 R10 尾界开始，自洽收紧 owner 半径至 1399/500，保留 1 列×3 行目录，共 234 原子、22 个 x、36 个 y。其完整 99,682 个候选盒枚举给宽度 2^−64 的真正最大端点包围，而非仅给一个下界见证。

另一个独立实现检查全部 168,498 个坐标 bbox（包含无法活跃的盒），得到

    C′ = X+sqrt(q),    ℓ=−sqrt(q),    q=319996899999/1280000000000.

9,410 个活跃盒中 9,393 个严格低于该端点，其余 17 个平手包含非紧 y bbox；不把平手次数说成互不相同的几何 owner。该确切根式位于小核 dyadic 区间内。完整轴证明见 `reviews/double-reflection-left-geometry-independent.md`。

## 3. 双缝闭合的全尺度尾界

完整实轴被长度 W_N 的相邻块分割，块内交替放原 N 列与其镜像；每个原胞元标签每块准确出现 N 次。由于 C>0、ℓ<0，N≥7 给

    W_N>6Px>30,    N<W_N/Px+1.

对长度 L=2r 的任意横区间：若 L<W_N，它至多穿过一个接缝，每标签点数至多 L/Px+2；若 L≥W_N，分成 q 个完整块和至多两个端块，点数至多

    (L−qW_N)/Px+2+qN ≤ L/Px+L/W_N+2.

纵向每标签点数至多 2r/Py+1。令 ρ=m/(PxPy)，便有

    Π_N(Q_r)/(4r²) ≤ ρ(1+Px/r)(1+Py/(2r)),           r<W_N/2,
    Π_N(Q_r)/(4r²) ≤ ρ(1+Px/W_N+Px/r)(1+Py/(2r)),    r≥W_N/2.

在 15≤r<W_N/2，第一式由 r=15 的准确值 0.99508438258…<1 排除。在 r≥W_N/2，第二式至多

    ρ(1+3Px/W_N)(1+Py/W_N)
      < (3ρ/2)(1+Py/(6Px)) = 0.9993268958951075… <1.

两个严格不等式均实际用 Fraction 核验。由此全部 r≥15 统一被排除，包含跨任意多个横周期的长半径。没有用两套完整原周期的错误密度 2ρ 支配，也没有漏掉每块完整 N 列与 W_N/Px 的差值。

## 4. 完整边界状态与系数恒等式

原条带由完整原周期支配，所以 r≥10 不活跃。又有 2Px>10。R46 的左右边态面积 E_L、E_R 分别使用完整窗口 [−2Px,2Px]×[0,Py] 与右端相同相位的平移窗口。这两个边带以及中间 N−4 个周期胞元完整覆盖原条带 halo，包括源块外侧接收区。因此

    A_strip(N)=(N−4)A+E_L+E_R=NA+b_x,
    b_x=E_L+E_R−4A.                                          (2)

双缝中的每次短半径 r<15 至多跨一缝，因为 W_N>30。原条带与右镜像的 halo 各自落在相邻长度 W_N 的闭块内；除零面积边界外不交。左邻镜像满足 τ_(−2W_N)T_(r_N)=T_ℓ，故其左缝是指定的真正左轴，不是另拼一个未对齐相位。

新增区域在每条缝两侧距离 <15 的带内。由轴反射对称，各完整接缝新增面积为半窗增量的两倍。左半窗为 [ℓ,ℓ+15]×[0,Py]，右半窗为 [r_N−15,r_N]×[0,Py]。原子距离接收中心 <15，所以每半窗只需接缝深度 <30 的源目录。N≥7 时，左目录 x<ℓ+30<30<NPx，不见省略的 i≥N 列；右目录 x>r_N−30>6Px−30>0，不见 i<0 的源。两半平面状态准确独立于 N。

所以每个双胞元完整面积为

    2(NA+b_x)+2γ_L+2γ_R.

除以 2Nm 即得 (1) 的恒等式。这里没有上下角，也没有 1/N² 项。γ_L 的 old 部分必须重新在 [ℓ,ℓ+15] 裁切，不能拿 E_L 代替。

该边界态恒等式不依赖 R41 的特殊质量：对任意具有相同尾界、局部性及分离条件的有限胞元周期输入，完整的 A、E_L、E_R、γ_L、γ_R 都给同一公式。

## 5. 新左缝完整计算与真轴误差

在 C′ 的有理外移轴上实际分别计算 old/new 全部接收 bbox：

| 状态 | 源原子 | 源坐标数 | 完整候选 bbox | 完整面积近似显示 |
| --- | ---: | --- | ---: | ---: |
| leftseam-old | 3625 | 209×62 | 39,259,817 | 79.01825718356… |
| leftseam-new | 8840 | 374×85 | 220,824,236 | 80.05106286250… |

两任务使用 40-bit dyadic 内外根界、准确前缀质量和整数矩形并面积，全部实际退出码 0。old 约 14.58 秒，new 约 77.06 秒；Python 3.13.18，macOS arm64，单进程，无随机数，无浮点容差。源代码使用冻结 R46 文件字节；本轮未改旧核。

将轴居中后，窗口固定为 [−15,0]×[0,Py]，原/镜像 x 端点各变化至多 δ=2^−64。任意 bbox 的接收矩形面积对称差至多 2δPy，包括原本为空、刚激活及窗口裁切的情形。共同目录采用轴区间较低端点，是整个区间的源标签超集；源顺序不变，因为该端点严格超过所有原源 x。

old/new 的全部 bbox 数上界分别为 42,858,585 与 256,306,875（不只数当前活跃盒）。两者共同给

    |Δγ_L| ≤ 2626724476943898165782541/14757395258967641292800000000000000
            = 1.7799377402646403…×10^−10.

已有 γ_R 同样采用 R46 真轴误差。主账本 `double_reflection_coefficients.py` 与独立 Fraction 账本逐分数一致，均严格证明 K<−12/5。没有用少数混合见证代替新 halo，也没有对多个 N 拟合或扫描。

## 6. 有限 L¹ 与正增益接口

完整证明见 `proofs/double-reflection-finite-l1-transfer-independent.md`。固定 N 后，Π_N 的胞元宽 T=2W_N>60，质量 M=2Nm。令 g_ε 是每原子由质量守恒小方块替换所得周期密度；均匀尾界 c<1 给 r≥15 时

    average_Q(g_ε) ≤ c(1+ε/15)²<1

（本轮 ε≤1/1024 足够且已准确核验）。有限局部 bbox 的内缩/外扩接收矩形给每胞元 halo 面积 S_ε→S。将源截成 Kx×Ky 个周期胞元再平滑，保留外侧完整 halo，并用横向 2、纵向 4 个 guard 胞元，得到上下因子

    (Kx−4)(Ky−8)/(KxKy),    (Kx+4)(Ky+8)/(KxKy).

取 Kx=Ky=j+9、ε_j=2^(−j−10) 后，有限 L¹ 输入 f_j 满足 J(f_j;1)→S/M=J(Π_N;1)。于是本轮负差距至少一半最终保留。一般推论是：只要同样完整边界账本给 K>0，就能通过该序列严格提高一般 C₂ 下界，增量可任意逼近 K/(mN)；这才是反射闭合的合法构造接口。

还得到一个一般近极值必要条件。若某周期模板的完整比值 J_per≥C₂−ε，且其双缝闭合在 N=N₀ 合法，则有限 L¹ 传递与恒等式给

    (K)_+/(mN₀) ≤ ε.                                         (3)

这是可计算边界缺陷对近极值的必要限制，保留全部重叠源集合和共同半径。它不能反向推出足够条件，也不证明从任意起点做反射迭代会趋于 C₂。后续递归机制必须改变并实际结算新的完整边界状态；继续增大同一 R41 的 N 无法改变已证负系数。

## 7. 证据与未解问题

主算轴、old、new 三项由 operations 顺序实际完成。几何审者另写独立左轴全 bbox 比较与共同目录计数；极限审者独立聚合准确系数并证明全尺度 finite-L¹ 传递。来源索引与白名单 manifest 记录实际源码、输入、输出 SHA。对冻结 R46 大状态的读取复用不计作本轮重跑；本轮也没有第二套全大面积核重跑。

本轮已决定规范 R41 双反射闭合路线的符号；后续应利用 (3) 组织新的近极值边界预算或其他递归构造，不再重复同模板 N 搜索。一般 C₂ 精确值、匹配的一般上界和构造性最优极值列仍开放。正式归档提交及远端核验由 root 执行，本工作区不提交或推送。


## 8. 最小复现命令

从整合后的课题目录执行；原 R41 输入、R46 大状态/系数和冻结算核作为已有依赖，R54 新输入由白名单提供。下列是同一输入与参数的复现命令，不把后来的复现声称为已经运行。本轮实际三个大状态由 operations 使用 Python 3.13.18 运行，其命令和退出码另存执行清单；核心聚合及误差界实际使用 Python 3.15.0。

```sh
python3 -B computation/reflection_edge_seam.py computation/inputs/r054_r041_left_reflected.json --output computation/results/double-reflection/r041-left-axis.json
python3 -B computation/reflection_edge_states.py computation/inputs/r054_r041_left_reflected.json --state seam-old --seam computation/results/double-reflection/r041-left-axis.json --output computation/results/double-reflection/r041-leftseam-old.json
python3 -B computation/reflection_edge_states.py computation/inputs/r054_r041_left_reflected.json --state seam-new --seam computation/results/double-reflection/r041-left-axis.json --output computation/results/double-reflection/r041-leftseam-new.json
python3 -B computation/reflection_axis_perturbation_bound.py computation/inputs/r054_r041_left_reflected.json --seam computation/results/double-reflection/r041-left-axis.json --output computation/results/double-reflection/r041-left-axis-error.json
python3 -B computation/double_reflection_coefficients.py --original computation/inputs/r041_final_0036.json --reflected computation/inputs/r054_r041_left_reflected.json --prior computation/results/reflection-edge/r041-coefficients.json --axis computation/results/double-reflection/r041-left-axis.json --old computation/results/double-reflection/r041-leftseam-old.json --new computation/results/double-reflection/r041-leftseam-new.json --axis-error computation/results/double-reflection/r041-left-axis-error.json --output computation/results/double-reflection/r041-double-coefficients.json
```

输出内的 elapsed、软件环境和绝对来源路径因运行环境变化，文件哈希可改变；应核对准确面积分数、源码/输入 SHA、参数与逻辑断言，不能把运行元数据不同误认为数学结果冲突。不可覆盖归档证书时，将以上输出改到 `/Volumes/data/CodexTemp/T2` 下的新运行目录，并同步修改后继输入路径。
