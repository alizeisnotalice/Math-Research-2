# R46 独审：边/角状态、反射轴扰动与全 N≥6 排除

日期：2026-11-15。独审工作区 `pcw-sol-reflection-limits`。本稿独立审读新 kernel、实际运行 17 个准确 toy states，并从 ops 的真实 JSON 用 Fraction 独立重算 R35/R41 的全部系数及轴扰动界；没有重跑这两个 204 原子模板的完整 halo 面积。

**结论范围：通过。** 对指定固定相位 R35/R41 模板，每个完整 N×N 胞元截断在 N≥6 后，关于真正最右 halo 边界反射所得比值都严格低于其自身真实周期 bulk 比值。这里还覆盖作者实际计算的有理向外轴。N=1,…,5、其他截断相位、其他模板或其他变换均未排除。不能据此声称 C₂ 上界、排除所有反射方法或排除某个未定义的更一般输入族。

主稿最终版本与新增有限 L¹ 误差已实际审读并绑定于文末；本稿同时绑定实际 kernel 和全部 state JSON 哈希。

## 1. 计算对象、裁切及相位逐项核对

审读的 edge kernel 为 `pcw-astra-reflection/.../computation/reflection_edge_states.py`，SHA-256 `ea664c7f9520b9983b13383fd032ff3110a6379b7dea381108ed49c88db9fede`；corner kernel 为 `reflection_edge_corners.py`，SHA-256 `53fe05c311c6601d4a7221e7518e3185b0b966ed737cd2a11c9916e5ce375fc3`。

本轮直接从代码核对以下对象，没有用文件 PASS 文字代替论证：

- bulk 的接收窗为 [0,P_x]×[0,P_y]，源为完整周期。left/right 接收窗为 [−2P_x,2P_x]×[0,P_y]，分别保留源列 i≥0/i≤−1；这正是左右各 4 个接收胞元，右态是将右边界按 NP_x 平移到 0，不能改成 i≤0。bottom/top 完全类似，源行 j≥0/j≤−1。
- 四个 old corner 的接收窗为 [−2P_x,2P_x]×[−2P_y,2P_y]，源象限分别是 (i≥0 或 i≤−1)×(j≥0 或 j≤−1)。每个角包括名义源窗口外的接收区域；不会把外侧角面积漏掉。
- seam 原源列为 i≤0，包含最后一个完整胞元；真正稳定轴的胞元坐标为 c_*，作者使用 c≥c_* 的有理上端。此处 i≤0 与 old right 边态的 i≤−1 是不同的参考平移，代码用法正确，不能强行统一。
- seam-left 接收窗为 [c−15,c]×[0,P_y]。因为 c≥真正 old halo 最右端，镜像旧 halo 完全在 x>c，故此窗中新 halo 面积减原 old halo 面积，准确等于左半纯混合新增 γ。全部纯混合中心都在 |x−c|<15，由对称性另一半同面积；不存在遗漏外侧接收区域或漏乘/多乘 2。
- seam-bottom/top 接收窗为 [c−15,c]×[−3P_y,3P_y]，原源 i≤0 且 j≥0/j≤−1，镜像保持纵坐标。这是上下各 6 个接收胞元的端状态。N=6 时底部行 −3,…,2 与顶部平移后行 3,…,8 刚好相邻，不重叠；N>6 中间有 N−6 个 bulk seam 行。

源窗按接收窗外扩所用 cap（old 10/new 15）生成，镜像用横坐标预像生成，没有只反射已裁好的旧 halo。二维 source bbox 的质量来自完整前缀和，枚举所有 bbox，检验共同半径 span 严格条件；不是把两轴半径独立优化。开边界在面积上可忽略，span 恰等于 2r 的框按严格阈值排除。

平方根 `rr=floor(D*r)` 用整数不等式复核；inner 采用 rr−1、outer 采用 rr+1，均在 cap 处截断。inner 比最紧 floor 更保守一格，不破坏有效性。receiver 与完整接收窗裁切后取并。未发现裁切、源相位、双倍或共同半径方面的错误。

## 2. 多项式系数与符号：无双倍遗漏

设 A 为完整周期胞元 halo 面积，E 为 left+right+bottom+top 四个边态面积之和，C 为四 old corner 面积之和。R10 的接收胞元 padding 为 2，因此

\[
A_N=(N-4)^2A+(N-4)E+C.
\tag{1}
\]

设 γ 为一纵周期左半 seam 新减旧面积，Γ 为底/顶两个左半端状态的新减旧之和。R15 纵 padding 为 3；轴已稳定且接缝带所需横深度已包含时，N≥6 给出

\[
D_N=2(N-6)\gamma+2\Gamma.
\tag{2}
\]

最右轴的更新独立证明在较小 N 已稳定；本稿的 N≥6 结论即使只用先前保守的 N≥5 轴稳定界也足够，无需把 N=2 的强化隐藏成前提。

令每胞质量为 m，定义

\[
L=E-8A+\gamma,
\qquad F_6=2E-32A+C+\Gamma.
\tag{3}
\]

有准确恒等式

\[
mN^2\left\{J(\mu_N+\mu_N^*)-{A\over m}\right\}
=L(N-6)+F_6.
\tag{4}
\]

等价地，常数项为 16A−4E+C−6γ+Γ。式 (3) 消掉了 N=6 时相同 γ 的相关误差，比把两个包含 γ 的区间分别组合更紧。只要严格认证 L<0 且 F₆<0，便直接排除所有 N≥6；不是只排除“足够大 N”。若只有 L<0 而没有 F₆ 符号，结论只能是充分大 N 不改善。

## 3. 实际 17 个 toy states 检查

`computation/review_reflection_edge_states.py` 导入两作者 kernel，但只在胞元 (10,10) 内一个质量 1 原子的极小准确输入上执行，不碰 ops 的实际模板计算。共调用 17 个 states，最大 25 个源原子、196 个 bbox；Python 3.15.0、单进程、Fraction、kernel bits=40，参考 √2 用 80-bit 定向有理根界。exit_code=0，结果 PASS。

独立解析预期为

\[
A=1,\quad(E_l,E_r,E_b,E_t)=(5/2,3/2,5/2,3/2),\quad E=8,
\]

\[
(C_{lb},C_{lt},C_{rb},C_{rt})=(25/4,15/4,15/4,9/4),\quad C=16,
\]

\[
\gamma=3/2-\sqrt2=:g,\qquad\Gamma=6g.
\]

实际 inner/outer 区间逐一包含这些准确值。式 (1) 还原 A_N=N²；式 (2) 还原 D_N=2Ng，故归一反射比值为 1+g/N，确认全部双倍和端点计数。

另把轴从 1/2 向外移 δ=1/100，严格解析值变为 γ=g+(2−√2)δ，实际 kernel 区间包含它，并与原 γ 区间严格分离。因此“轴误差很小所以可完全忽略”不是正确的证明；必须实际控制其影响。实际结果为 `computation/results/reflection-edge-state-review/exact-toy.json`。

## 4. 真正轴与有理向外轴：完整共同 catalog 的误差证明

作者轴证书给 c_*∈[c_lo,c_hi]、有理轴 c=c_hi、0<c_hi−c_lo≤2^−64，且 c_lo>全部原胞元横坐标。独审实际核对了这些字段和所有 state 的轴/输入/seam 哈希一致性。

不能只计算新 mixed 矩形的物理端点移动而漏掉裁切窗 [c−15,c] 的移动及 old slice 的变化。有效方法是使用轴居中坐标 t=x−c，将接收横窗固定为 [−15,0]。原源坐标 α+iP_x−c、镜像坐标 c−α−iP_x 随轴移动至多 δ，且 source x 排序不变。每个完整 bbox 的质量、根半径及 y 端点不变，两个 x 端点各移动至多 δ。故该接收矩形的面积对称差至多

\[
2\delta\min\{H,2R\},
\tag{5}
\]

H 为该接收窗高度。包含空框激活、横裁切、端点相等等情形，不要求 chamber 内部。

lower-axis 共同 catalog 确实覆盖整个轴区间：原源只可能由条件 x≥c−15−R 进入窗；镜像的原坐标预像只可能由 x≥c−R 进入窗。另一侧上界因 x≤x_max<c_lo 恒不活跃。轴增大只会删源，不会加入 lower-axis catalog 之外的源；纵向范围与半平面限制不随轴改变。保留 lower-axis 全部源并随轴移动，仍可在目标窗内准确表示真实 halo。

用该共同 catalog 的全部 x/y bbox 数

\[
K={n_x(n_x+1)\over2}{n_y(n_y+1)\over2}
\]

乘 (5)，不会遗漏当前 inner/outer 中不存在但在轴区间内出现的框。对 new 与 old 两者都加误差；再分别对底顶端相加。不能拿 surviving unique rectangle 数直接代替 K。

作者 `reflection_axis_perturbation_bound.py`（SHA-256 `63b932fbe200af7d301657c05906d23196d3d6cc1cd5ee6f72a5ace8a6402bad`）采用的正是这个方法。独审另写一种源生成实现，只用上述必要的原坐标 lower bounds，并独立对 lower/upper 两端标签集合检查单调包含，实际复算全部 nx、ny、K、高度和 Fraction 误差。R41 六态全部逐项准确相等，得到

\[
|\Delta\gamma|\le1.883153207\times10^{-10},
\qquad |\Delta\Gamma|\le1.582171216\times10^{-9},
\tag{6}
\]

这里显示的是向上放宽十进制，准确分数见独审 JSON。R35 由独立 catalog 直接计算，误差约为 1.8258210×10^−10 与 1.5339394×10^−9。

**第二条不依赖实际 catalog 的粗证书。** 每标签原 x 共同范围宽小于 30，可保守取 6 列；bulk y 可保守取 8 行，end y 可保守取 13 行。204 标签给出 nx_old≤1224、nx_new≤2448、ny_bulk≤1632、ny_end≤2652。把这些纯标签上界代入全部 bbox 计数和 (5)，得到共同适用于两模板的粗误差

\[
|\Delta\gamma|<10^{-5},\qquad |\Delta\Gamma|<10^{-3}.
\tag{7}
\]

独审实际 Fraction 值约为 2.970881886×10^−6、9.411764731×10^−5；即使只用这条明显更宽的独立粗界，最终符号仍不变。

## 5. ops 真实状态独立聚合与严格符号

真实 state JSON 位于 `pcw-operations-terra/reflection/r0XX-*.json`；R35 的 bulk 取冻结 `pcw-vertical-compute/.../results/round035/primary_old_R10.json`。每模板共读取 15 个面积证书，逐一检查输入 SHA、代码 SHA、union-module SHA、质量、阈值、轴/seam SHA、inner≤outer 及 tail<1，然后独立做区间运算。所有文件哈希清单写入对应独审 JSON。

`computation/review_reflection_r041_coefficients.py` 实际重算 R41 全部系数，未信任作者汇总的小数。`computation/review_reflection_r035_coefficients.py` 对 R35 完成同一独立工作；R35 还独立生成自己的 axis catalog，不假设其位移对 K 没影响。两个进程均 exit_code=0、PASS。

加上真正轴的扰动误差后，区间上端近似为：

| 模板 | L 真轴区间上端约值 | F₆ 真轴区间上端约值 |
| --- | ---: | ---: |
| R35 | −9.105721697777488 | −53.293751299325 |
| R41 | −9.105721600526445 | −53.2937507507865 |

严格认证实际使用分数而非这些近似显示。两模板均验证干净的有理结论

\[
L<-9,\qquad F_6<-53.
\tag{8}
\]

并且用第 (7) 条纯标签粗界重新扩大误差后，(8) 仍通过。于是对每个 N≥6，由 (4)

\[
mN^2\{J_{\rm reflected}-A/m\}
<-9(N-6)-53<0.
\tag{9}
\]

它排除的是超过该指定输入**真实周期 bulk**比值；不应将真实 bulk 与当前有限精度打印的某个有理下界混为一谈。没有计算 N≤5 的角重叠过渡态，也没有给出其他相位的排除。

## 6. 实际核验范围、文件与复现

本轮实际数学运行全部为 Python 3.15.0 `/Users/zzh/.local/bin/python3`，单进程、Fraction/整数根、无随机种子、无浮点容差。17 toy states 是本轮实际 area kernel 调用；R35/R41 的大状态面积是读取 ops 证书后独立聚合，不能写成本轮重跑。共同 catalog 只生成/计数，未枚举其大 bbox 面积。没有 MATLAB、大搜索或额外占用 ops 作业。

必要附件：

- `proofs/finite-state-reflection-audit.md` 及其 `computation/finite_state_reflection_audit.py`、`computation/results/finite-state-reflection/exact-check.json`：有限状态阈值与小反例。
- 本独审 `reviews/reflection-edge-state-review.md`。
- `computation/review_reflection_edge_states.py`、`computation/results/reflection-edge-state-review/exact-toy.json`：17 个真实 toy state 调用。
- `computation/review_reflection_r035_coefficients.py`、`computation/results/reflection-edge-state-review/r035-independent-coefficients.json`：R35 实际证书绑定及独立 Fraction 系数/轴误差。
- `computation/review_reflection_r041_coefficients.py`、`computation/results/reflection-edge-state-review/r041-independent-coefficients.json`：R41 同上，并逐项核对作者误差证书。

三个 review 脚本的默认源路径明确绑定本机这些工作区；复现时应保留其 source 目录，或据归档来源索引重建相同路径。输入 JSON、state JSON 和作者 kernel 的确切 SHA 都保存在结果文件内。按各脚本绝对路径执行 Python 即可复现本轮审计；不自动启动 ops 原大面积计算。

冻结 R30/R37 附件未改。本轮新增材料未提交/推送，交 root 白名单整合。一般 C₂ 的准确值、匹配一般界、一般 L¹ 达到性仍开放。

## 7. 主稿最终版本绑定与有限 L¹ 半亏损

已实际全文审读主稿 `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-reflection/planar-centered-square-weak11/proofs/high-ratio-reflection-edge-obstruction.md`，大小 13870 bytes，SHA-256 `9bf263dc2914990cd2643dd35f043aa29c5ac83f5ee817257f3b267c7ed9f8c8`。其固定相位、真轴区间、N≥6 多项式、两模板符号和局限与本独审一致，最终主稿通过。本次没有把 singleton 下界见证当成准确全局最大性的证明；主稿 §3 已明确使用 geometry Sol 的另一个 63,000 bbox 严格比较证书，本独审的 L/F₆ 结论本来就对完整轴区间稳健，不依赖准确 singleton 等式。

主稿新增 §8 的有限 L¹ 常数正确。反射后原子质量 T=2mN²、原子数 n≤408N²，令 K=2ⁿ−1。每个原子质量支撑于其中心半边长 η 方块时，严格 halo 包含于全部原子子集接收矩形向外扩 η 的并，故完整面积增加至多

\[
E_\eta\le K(4\eta\sqrt T+4\eta^2).
\]

原来为空的侧长使用正部后仍满足同一面积增量界，因此没有漏掉激活框。若 η≤1 且

\[
\eta\le{9N-1\over4K(\sqrt T+1)},
\]

则 E_η≤4Kη(√T+1)≤9N−1。原原子比值亏损是 (9N−1)/(mN²)，对应**面积**亏损 T×该值=2(9N−1)。所以平滑增加面积至多 9N−1 后仍严格低于 bulk 至少一半原比值亏损；严格性来自原子不等式本身。不存在漏掉因子 2。

此 L¹ 结论准确涵盖“对反射后的原子输入逐点使用相同 η 方块平滑”，也等价于先同尺度平滑原子两半并保留原原子反射轴；它没有重新计算平滑后原半输入的真正 halo 边界，也未排除任意其他平滑尺度或重新选轴的 L¹ 族。主稿未将其扩展为后一种断言。

另外实际运行 `computation/review_reflection_l1_half_gap.py`，对 N=6、7、10 用 √T 的严格有理上界选择正有理 η，准确验证 E_η≤9N−1 及至少半亏损，exit_code=0，PASS。该核验只计算整数/Fraction 常数，最大 n 上界 40800，没有生成任何原子或平滑大文件；一般参数结论由上面的解析证明承担。输出为 `computation/results/reflection-edge-state-review/l1-half-gap.json`。

本独审、全部脚本与实际结果现已冻结；最终归档源/目标/大小/哈希见同目录 `r046-sol-reflection-limits-manifest.json`。冻结 R30/R37 六文件未改；正式提交及远端推送由 root 执行。

## 8. 切相位补充的最终绑定

已实际审读新补充 `proofs/high-ratio-reflection-phase-supplement.md`，2133 bytes，SHA-256 `8730f0c3a6ce644b1ba0ad8262a8cde6cfcc3c88571ea884faf567f3fdc91e98`，原冻结主稿 `9bf263dc…` 保持原字节。相位补充通过，不改变任何大面积状态、系数或 N≥6 结论。

其充分条件需连同真实相位写为 C+4P_x>30，而不是只引用 5P_x>30。独审在本轮另用本机 Fraction 对两份冻结输入实际计算，均准确得到

\[
x_{\max}+4P_x-30
={424640887836625294630001\over150000000000000000000000}>0.
\]

由 C>x_max 得 c_N−30=(N−5)P_x+(C+4P_x−30)>(N−5)P_x。所有编号≤N−6 的源都严格位于 (N−5)P_x 左侧，故不能进入接缝窗口任何 r<15 的严格见证；N≥5 时同样看不到 x<0 的左截断。镜像和半径更小的 old 状态均被涵盖。纵向 padding 和 N≥6 端状态分解因此完整成立。

同时实际审读并绑定 `computation/reflection_cut_phase_check.py`，2032 bytes，SHA-256 `262cd7676c510d26c348bfb33e60d98c387efe8090cd2a7db57a9001b4a79bd7`；其实际结果 `computation/results/reflection-edge/cut-phase.json`，1520 bytes，SHA-256 `2f2a1313d889ebfb28b825a8607eba2d75137f4ef379cd76f984d8ccd84af152`。作者检查的输入哈希、半开横相位与准确余量均一致。独审未重跑任何大状态；本轮新增实际计算只有上述两份输入的 Fraction 余量复核，全部 PASS。

最终 manifest 现列 12 个独审附件及这 3 个补充附件的绝对 source/dest、SHA-256 和字节数；所有其他冻结 payload、主稿以及 R30/R37 附件保持原字节。
