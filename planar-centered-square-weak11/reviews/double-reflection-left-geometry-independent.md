# R54 真左轴、左缝相位及双闭合系数独立审查

2026-11-15；sol_reflection_geometry。仅写本人工作区；R30/R37/R46 冻结文件未改，无提交推送。审查对象为冻结 R41 canonical cell，二维轴平行中心共同半径正方形，严格阈值 1，全尺度。沿 y 保持 Py 周期，沿 x 取 N 列后在真正左右 halo 轴双反射闭合。本文新增完整左轴单根式与左轴扰动目录独审；未重跑大型 seam 面积。左右两侧半缝增量须区分，不能套用四分之一角增量。

## 1. 输入与安全的取反相位

原 R41 输入 SHA256 `ab714c7f457ccf69f68fd629c947021faab16631c8995caa16b67c1ca2752db3`；变换输入 `r054_r041_left_reflected.json` SHA256 `deb65c6011d01f61253cb1db6b4b0e919f014ad32711677b9ee3fd80645694bf`。实际逐项读取 204 点并用 Fraction 验证

X=max(alpha)=979821360288625294630001/150000000000000000000000，min(alpha)=0，X<Px，alpha'=X-alpha，beta'=beta，w'=w。

所有 alpha' 已在 [0,Px)，没有取模重排。故 X-(alpha+iPx)=alpha'+(-i)Px，原左半平面 i>=0 精确映为变换右半平面 i'<=0，原 N 列映为 1-N<=i'<=0。这与先取 -alpha 再各标签取模造成不同 cutoff 的错误不同。任一 centered square 经 x'=X-x 仍是同一半径 square；阈值、质量、面积和 Py 不变。

若 C' 是变换右半平面的真右 halo 轴，则原左半平面的真左轴 ell=X-C'。原窗 [ell,ell+15]×[0,Py] 精确映为 [C'-15,C']×[0,Py]，old/new 两个面积各自保持，gamma_L 不多乘 2。关于 ell 的反射也共轭为关于 C' 的反射。operations 使用 Cbar'>=C' 时对应原轴 ellbar=X-Cbar'<=ell，方向确为向外。

## 2. 新左轴的完整目录和独立单根式证书

本人新脚本 `computation/double_reflection_left_axis_independent.py` 只读原输入、变换输入及 operations 新轴输出 `r054/r041-left-axis.json`。没有导入其算核，没有使用其 bbox survivor/根下界剪枝。先逐轮重新按输入重建 source 标签行，核对全部 bootstrap 质量 B 和严格正差 4Rnext²-B。各轮从 R10 收紧至 1399/500，最后目录是单横列、三纵行、234 原子，22 个 x 和 36 个 y。目录完整性来自右端竞争 u>X 蕴含最左 source >X-R、横列 ceil(R/Px)、纵 row ceil(2R/Py)+1；无限纵周期允许将最低 row 平移到 0。不是把有限纵行当总输入。

对完整 168498 个坐标 bbox 作准确 prefix 质量及 span²<M 比较：5775 零质量、153313 空接收、9410 非空。取 m0=319996899999/320000000000，q=m0/4。每个非空 bbox 最左横坐标 a、p=M_B/4、g=X-a 满足 a+sqrt(p)<=X+sqrt(q)：839 个 p<=q，8571 个 p>q 但 p-q-g²<=0，无正 L 平方分支；9393 严格，17 取等。17 个取等 bbox 的横坐标均为 X，质量均为 m0；部分 y 边界不紧，不能要求它们全是几何 singleton bbox。其实际 source singleton 当然实现该最大端点。

因此严格得到

C'=X+sqrt(q)，ell=-sqrt(q)，q=319996899999/1280000000000。

q 非有理平方。另准确平方核验 operations 的 [lo,hi] 包含此真根，hi-lo=2^-64。既有旧 144 原子证书不适用于这份新左输入，本次明确重新完整检查 234 目录。输出 `results/double-reflection-left-independent/exact-axis.json` 保存输入/脚本 hash、全部计数和严格余量。Python 3.15.0，Fraction/isqrt，单进程，无随机数/浮点容差，最终 exit_code=0；首次运行只因过严的 y bbox tie 断言失败，放宽到实际需要的相同横端/质量后全部证明断言通过，没有改变比较判据。

## 3. N>=7 的有限列稳定与局部源完整性

N 列纵周期条带 sigma_N 的旧 halo 左轴就是 ell，右轴 r_N=(N-1)Px+C，其中 C 是 R46 已证原输入右半平面真右轴。新左目录只用横列 0，故每个 N>=1 已包含其轴 owner，反向由子测度支配给等号；有限纵行并非限制，因为 sigma_N 纵向本来无限周期。

严格有 ell<0、C>X>0。N>=7 给 W=r_N-ell>(N-1)Px>=6Px>30。左半缝窗内任何 r<15 square 只见 x<ell+30<30，而被 N 列截去的原右侧列 source 均 x>=NPx>=7Px>30。镜像预像相同，故左缝精确等于无限左半平面状态。

右半缝窗内任何 r<15 square 只见原侧 x>r_N-30>(N-1)Px-30>=6Px-30>0，被右半平面补入的负编号原 source 均 x<0，故有限条带右状态同 R46。两方向所需的严格有理余量 6Px-30、7Px-30 均实际 Fraction 检查并记录。W>30 同时排除一个 r<15 square 穿过两条相邻缝；两个宽15新增带不交。因此左右 gamma 都不依赖 N。

## 4. 低轴共同目录及真轴转移

审读冻结 `reflection_axis_perturbation_bound.py` 后，本人另写 `double_reflection_left_perturbation_independent.py`，独立计算 inclusive source lattice ranges，不导入冻结 source_atoms。固定轴居中窗 [-15,0]，原 source x-c、镜像 source c-x；lo>X 保证两侧符号与全排序固定，移动 c∈[lo,hi] 的每个横坐标至多 delta=2^-64。

原 source 的预像横窗为 [c-15-R,c+R]；镜像为 [c-R,c+15+R]，且所有原 source <=X<lo。故真正影响裁切窗的下 cut 随 c 增大收紧，lo 目录是整个轴区间的 superset。实际核对 hi 的标签目录包含于 lo 目录。共同目录允许保留后来离开的 source；它们离接收窗大于 R，不会生成非空 cap-R 接收矩形。

独立结果：old R10 目录 3625 原子、209×62 坐标、全部 Kold=42858585 盒；new R15 目录 8840 原子、374×85 坐标、Knew=256306875。它们是全坐标盒计数，不是 active survivor 数。每个 bbox 的质量和 cap-root 不变，横端点移动 <=delta，裁切也 1-Lipschitz。包括空盒、激活、裁切墙，区间对称差长度 <=2delta；高度不变且 <=Py，故面积变化 <=2delta Py K。并的对称差由各盒之和控制。

两态合计准确得到

|gamma_L(true)-gamma_L(outward)| <=2626724476943898165782541/14757395258967641292800000000000000，

显示约 1.7799377402646403e-10，与作者新轴误差记录准确相同。原坐标轴移动方向相反不改变绝对误差，无额外 2。本文只认证 seam old/new；纵周期 R54 不需要 corners，未独立重算作者附带角误差。独立输出 `results/double-reflection-left-independent/perturbation.json`，实际 exit_code=0。

## 5. 全尺度尾界及准确边缘公式

既有 `double-reflection-closure-audit.md` 的全尺度分块证明通过：横 period 2W、纵 period Py、每胞元质量 2Nm。每个 W 块含 N 个该标签源；W>(N-1)Px 给 N<W/Px+1。长度 L<W 区间只跨至多一缝，点数 <=L/Px+2；L>=W 时分完整块和两个碎块得 <=L/Px+L/W+2。不能把两份原周期密度相加为 2rho。

因此短半径用 rho(1+Px/r)(1+Py/(2r))，长半径 r>=W/2 用 rho(1+3Px/W)(1+Py/W)<rho*(3/2)*(1+Py/(6Px))<1；r>=15 的短段也由 R46 R15 界排除。后一严格不等式本次另实际由原质量/周期核对；不依赖有限 N 扫描。故所有共同半径、所有中心与全部尺度都在 cap15 目录内，而非改了算子尺度。

旧纵周期条带面积每 Py 为 (N-4)A+E_L+E_R=NA+b_x，其中 b_x=E_L+E_R-4A；没有上下边/角。一个 2W 胞元包含两份旧条带 halo，其新增是两条完整缝，分别 2gamma_L 与 2gamma_R。因此

|H(Pi_N)|cell=2NA+2b_x+2gamma_L+2gamma_R，
J_N=A/m+(b_x+gamma_L+gamma_R)/(mN)，N>=7。

gamma 是半缝增量，两个反射半边在总面积中各乘2，除质量2Nm后恰抵消。这里不是四分之一 patch 增量，没有 1/N² 角项，也不除胞元几何面积。采用外移轴时周期必须同步取2(rbar-ellbar)；回真轴只需给两个 gamma 的扰动，不额外增加 period-area 项。

## 6. 审查状态与尚缺计算

新左输入相位、真左轴、N>=7 局部状态稳定、低轴共同目录、gamma 扰动界和不依赖 N 的公式均通过。冻结 R46 bulk/edge/右缝结果明确复用；本轮未重新计算其大型面积。新 leftseam-old/new 均已由 operations 实算完成，作者主账本与 limits 独立账本准确聚合得到真轴 gamma_L≈[1.0328056783975574,1.032805679474352]，K≈[-2.436023761291603,-2.4360237578629964]，严格 K<-12/5。小数仅为显示。该符号由主/limits 的准确聚合承担；本人本次读审其结果和几何接口，没有新增大型面积重跑，也不把读取记录记作本人重新执行系数脚本。故规范 R41 相位、阈值1、每个 N>=7 的双闭合族严格低于原周期 bulk 比值；一般 C2 仍开放。

归档白名单、来源 hash/大小及建议目的路径见同目录 `results/double-reflection-left-independent/manifest.json`。所有新文件均小于10MiB；main/T 保持只读，root统一整合。

## 7. 最终主证明 §§2–5 的一次定向读审

实际完整读取作者 `proofs/double-reflection-periodic-strip-obstruction.md`，本次定向审核 §§2–5。读取版本 SHA256 `6ef36428cc0e5fa68707327987b2d553185108b9b8ff8347b4d6e64eea685b09`，12440 bytes。同步读取主账本 `computation/results/double-reflection/r041-double-coefficients.json`，SHA256 `1e70141ac81fe9f71930e75e3c329686f3df6a0b8268c00dd41ae819f795711e`，6614 bytes；本次不重新执行其聚合脚本。

结论：通过。§2 不取模变换、234目录和两种 bbox 计数区别与本独审一致；§3 长短半径分块含多胞元长尺度且 leading 密度不翻倍；§4 N>=7 左右局部性、两份旧 halo 不交及完整缝各为半缝两倍的面积账本均正确；§5 低轴共同目录、空/激活 bbox 扰动和真轴转移与本人实际独立核验完全一致。old/new 40-bit 面积来源与较宽完整 bbox 误差计数明确区分，没有把 39,259,817/220,824,236 的预剪候选计数误作全部扰动目录盒数。

最终系数断言的范围明确为 R41 canonical相位、真左右轴、阈值1、N>=7；本文不增加所有阈值、R35或任意平滑族结论。主证明 §6 的 finite-L1 定理另由 limits 审者负责，本次定向几何审查不冒充该附件的额外独立重证。
