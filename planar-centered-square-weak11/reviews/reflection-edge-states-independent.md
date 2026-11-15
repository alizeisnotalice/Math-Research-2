# R46：局部边界/角核、真实轴及R41反射族的独立审查

日期：2026-11-15；sol_reflection_geometry。对象固定为R35/R41两个冻结周期输入、N×N cell截断、一次关于最右轴或其认证有理外移轴的横向反射、严格阈值1。结论：局部核覆盖与上下包围正确；真实最右轴有独立单根式证书；两模板全部profile的独立准确组合及轴误差证明，每个N≥6的反射比值均低于原准确周期bulk比值超过(9N−1)/(mN²)。详细局部核测试以R41为基准，双模板最后绑定及L¹外扩上界见§§9–10。N≤5、其他模板、其他反射策略、任意云团变体及一般C2均未被排除。R30/R37冻结材料未改。

## 1. 绑定输入、核与证据性质

输入R35 SHA256：5021c54704af45b3ef8289c67a245926473434640372a934b2f15c4485f0c026；输入R41 SHA256：ab714c7f457ccf69f68fd629c947021faab16631c8995caa16b67c1ca2752db3。完整局部核实际读取并做小例执行的版本：

- `reflection_edge_states.py`：6210 bytes，SHA256 ea664c7f9520b9983b13383fd032ff3110a6379b7dea381108ed49c88db9fede。
- `reflection_edge_corners.py`：6357 bytes，SHA256 53fe05c311c6601d4a7221e7518e3185b0b966ed737cd2a11c9916e5ce375fc3。
- 只导入其并面积函数的冻结`periodic_halo.py`：6119 bytes，SHA256 b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33。
- `reflection_axis_perturbation_bound.py`：SHA256 63b932fbe200af7d301657c05906d23196d3d6cc1cd5ee6f72a5ace8a6402bad。
- R41轴扰动实际输出SHA256：b6c56ea153aa8c2e8da8fc6611af62d7427f4a6aecfc56bc847ab6f2e7760230，3574 bytes。

本次完整审读核与一般分解证明；独立执行15种source-filter测试、15种可解析稀疏面积测试；独立写144源全bbox单根式检查；独立写15个实际面积JSON的Fraction组合检查。大型R41局部面积由operations实际执行，本次不将它们写成本人第二次大算。轴误差小核另实际重放，输出与作者JSON逐字节相同。上述区别均在实际记录中保留。

基础尾界/半径bootstrap/N9安全分区见本人新稿`proofs/reflection-periodic-tail-locality.md`。该稿先证明R15及144目录，本次进一步给精确C和更紧N6分解；不是把两套X/Z命名混用。

## 2. source_atoms和完整bbox算核

edge的left/right窗口为[-2Px,2Px]×[0,Py)，源指标分别i≥0、i≤−1。right应相对几何右边界NPx平移，因此i≤−1正确，不应与seam相对最后cell的i≤0混同。bottom/top窗口为[0,Px)×[-2Py,2Py]，源j≥0、j≤−1同理。

seam窗口为[cbar−15,cbar]×[0,Py)，源i≤0，镜像也是同一半平面关于cbar反射。镜像source-window预像[2cbar−(x1+R),2cbar−(x0−R)]准确；纵轨道没有反射。corner-lb/lt/rb/rt分别使用对应两个指标截断，窗口[-2Px,2Px]×[-2Py,2Py]。seam-bottom/top窗口[cbar−15,cbar]×[-3Py,3Py]，j分别≥0、≤−1，并对new状态加镜像。条件分支确实覆盖所有8角状态。

原source计数尾界R10严格<1；seam-new使用不翻倍leading密度的R15严格<1，旧seam仍用R10。source_atoms精确加入接收窗R邻域内所有闭边界源，正/负指标ceil、floor均正确。bbox中全部源质量由整数prefix求和，反射/原子偶有同位置也以质量相加而非丢失；本案两副本因轴在全部原source右侧实际无碰撞。

完整性：任一r<R见证捕获集合的bbox在该闭邻域内，遍历坐标端点及span≤2R覆盖它；bbox内其他原子也全部在闭source-window内，所以质量不漏。反向允许额外捕获，不要求bbox恰为唯一捕获集合。receiver只用t=min(R,√(M_B/4))；非空条件span²<M_B的准确平方判据正确，空激活边界不会被当正面积。

scale包含所有位置/窗口分母、D与R分母，weight_scale包含质量分母，整数转换准确。根界rr满足rr²B≤M A D²<(rr+1)²B，内根(rr−1)/D、外根(rr+1)/D方向正确，二者均cap到R；不能因为内矩形空就丢外矩形，代码分别保留。盒级x/y窗裁切与source span预剪只删除cap R后不可能进入窗口的盒。导入旧union_area处理的是矩形并，不重复加面积。故各state内/外分数分别是该完整局部halo面积的严格合法下/上包围。

## 3. 真正最右轴的独立单根式证书

先前半径bootstrap已把所有右端u>xmax候选归约到同一144原子有限目录（72 labels×两纵row），R35/R41在该目录完全相同。其20个x、24个y产生63000个闭坐标bbox。本次新checker没有使用operations的根界下界剪枝，也没有只枚举survivor或猜测owner。

取

xmax=979821360288625294630001/150000000000000000000000，

mright=124393118500658757956428526782554147336111/25000000000000000000000000000000000000000000，q=mright/4。

新checker独立取x=xmax处的最大singleton质量为mright。对任一非空bbox，令a为最小横端、p=M_B/4、d=xmax−a≥0。目标比较a+√p≤xmax+√q通常可按p≤q，或L=p−q−d²≤0，或L²≤4d²q准确判定。本案实际得到更强的全bbox有理不等式

M_B≤mright+4(xmax−a)²。

因为√p≤√(q+d²)≤√q+d，所有右端均≤xmax+√q。63000盒中3230零质量、54704接收为空、5066非空；其中64个p≤q，余5002个均L≤0，无需正L平方分支。5062个右端严格较小；只有4个x=xmax、同质量mright、y平移singleton取等号。它们确实是合法source，所以

C=xmax+√(mright/4)

是严格的真正最右轴，而非仅一个很窄区间。Fraction约分后的q并非有理平方（分子/分母整数平方测试），因此C无理。先前cN=(N−1)Px+C对N≥2成立，R35/R41共用此C。

`computation/reflection_axis_single_root_check.py`与`results/reflection-axis-single-root/exact-certificate.json`保存全部计数、等号见证和最小正有理质量余量。当前checker可仅凭相邻common144目录运行；`--seam`可选，用于额外检查operations的64位区间确实包含该真根。本次实际带--seam运行，exit_code=0，约0.067秒。没有大代数域或浮点符号判断。

## 4. 有理外移轴是合法构造

operations取cbar为已证明C区间的有理上端，宽度2^-64。因为C≤cbar，旧halo始终在x<C≤cbar，反射旧halo在x>2cbar−C≥cbar，两份旧halo仍不交。关于cbar反射是合法新输入；不能说cbar等于无理C。

对原source全部在cbar左、镜像全部在右的测度，半线横交长之和≤2r仍保证R15尾界，leading密度不翻倍。R41新旧state严格计算的是这个有理外移构造。回到真正C需要下一节的扰动证书，不能只说2^-64“很小”就丢掉差异。

## 5. N6分解、leading与角常数

旧R10下2Px>10、2Py>10。按左右各2cell的固定边带、上下各2cell的固定边带及其四角分割，所有N≥6有准确恒等式

A_N=(N−4)²A+(N−4)E+Ccorner，E=L+R+B+T。

这里A是完整周期cell halo面积，E是四个edge窗口面积之和，Ccorner是四个corner窗口面积之和。故

A_N=A N²+bN+d，b=E−8A，d=16A−4E+Ccorner。

新增面积只可能在最右轴左右距离15内：左半离轴≥15的中心，其r<15方形碰不到任何镜像source，new与old halo完全相同。纵方向取上下各3cell带，因为3Py>15；在N≥6时两带不互见，中央N−6个纵period完整。水平方向N≥6有(N−1)Px+C_lower>30，使该seam窗口的R15邻域不见原块最左截断或镜像最外端。故对任一c∈[C_lower,cbar]，包括真C与有理cbar，准确有

D_N=|H(νN)|−2A_N=2(N−6)γ+2Γ。

γ是同一seam窗口每纵周期的new−old面积；Γ为bottom/top两个seam角窗口new−old面积之和。top应按NPy平移后用j≤−1，角代码已如此实现。两项均≥0，但将式子展开成2Nγ+常数时常数未必非负。

不要把宽15的seam-old面积直接替代宽4Px的right edge面积：它们窗口不同。现行公式中b只由四个edge与bulk构成，γ只由同窗new−old构成，没有发生这种混用。

质量为2mN²，所以

Jnew−A/m=F_N/(mN²)，

F_N=(b+γ)N+K，K=16A−4E+Ccorner−6γ+Γ。

leading符号准确是b+γ，没有多一倍、漏除Py或重复除质量。更稳妥的全族检验为

F_N=F6+(N−6)(b+γ)，F6=2E−32A+Ccorner+Γ。

N=6处γ严格抵消，组合上下界时应先化简再代区间，避免重复的不必要放宽。负leading本身只能排大N；还要认证F6<0才能排全部N≥6。未算角时可先用Ccorner≤64PxPy、Γ≤180Py，得K≤80PxPy+180Py，但这不是已计算角系数。本次实际角输出齐备，已使用准确区间，见§7。

## 6. 真实轴扰动的共同目录界

用轴居中坐标固定接收窗[-15,0]。当c从lo变到hi，原source横坐标x−c向左移动，镜像source c−x向右移动，各≤δ=hi−lo。lo>xmax保证所有原source负、所有镜像source正，因此全部横坐标排序在整个轴区间固定，重复坐标也始终重复；纵坐标不动。

共同目录取下轴lo处接收窗的R邻域内所有source。对原source，进入该邻域的下cut为x≥c−15−R；对镜像source则为x≥c−R，二者都随c增大而收紧，所以lo目录确为整个区间的superset。不是取某个axis的survivor矩形列表当共同覆盖。

固定共同横/纵坐标索引bbox后，所含标签和质量不变，临界根与cap R不变；横接收端点至多移动δ，裁切到固定[-15,0]仍1-Lipschitz。纵高度不变且≤min(window height,2R)。即使跨激活/空盒/裁切墙，区间对称差长度仍≤2δ，因此每bbox面积对称差≤2δ min(height,2R)。矩形并的对称差由全部bbox之和支配。

小核使用完整K=nx(nx+1)ny(ny+1)/4，不用span/root surviving计数；包括old和new两边，六个状态全部付费。因此轴扰动输出严格给

|Δγ|≤约1.8831532061×10^-10，|ΔΓ|≤约1.5821712159×10^-9。

准确分数在作者`r041-axis-error.json`中。本次逐行重核K及2δ·height·K，另实际重放该小核约0.114秒，输出与作者3574 bytes JSON逐字节相同。新旧halfplane各自有统一R10/R15尾界，故共同cap公式对整个轴区间都合法。

## 7. R41全部N≥6的准确组合结论

本人独立读取operations的R41全部15个state JSON，核状态名、输入hash、source kernel hash、union_module hash、seam hash、tail及内外面积方向。每个实际文件的完整hash/体积保存在新的`computation/results/reflection-edge-kernel-independent/r041-assembly.json`，不是只读取浮点字段。

有理外轴时b约−10.1481833957、γ约1.0424617934，leading约−9.1057216023。独立Fraction先准确组合各区间，再加上述真实轴误差，得到真正C的合法上界

b+γ(C)<−9.1057216005，

F6(C)<−53.2937507507。

这些小数只作阅读显示；实际判负使用JSON中的准确负分数。因为F_N(C)=F6(C)+(N−6)[b+γ(C)]，对全部整数N≥6均严格<0。因此固定R41原子模板的N×N截断，作一次真正最右反射后，仍严格低于其准确旧全周期bulk比值A/m。有理外移构造也同样低于。

范围限定：这不是一般C2上界，没有核N=1,…,5，没有排除其他反射轴、纵向反射、多次反射或改过质量/坐标的模板，也不对任意非零云团尺度的L¹变体给上界。两输入虽共用C，本段全族数值组合只读取R41全部15state，本节R41数据不自动覆盖R35；R35已按§9另行准确复核。

## 8. 实际独立小例和归档建议

`reflection_edge_kernel_smoke_independent.py`实际对15种source_atoms状态，与独立有限指标穷扫的闭窗口字典逐项比较；再对15种state以周期4×5单单位原子的稀疏格作准确面积检查。有理外轴1时R1尾界7/8<1，任意两source的L∞距离≥2，故全部halo都仅由singleton形成，真面积无需待审bbox核即可解析给出。两核的dyadic内外区间逐一包含这些真面积，全部PASS，约0.040秒。该测试同时涵盖i≤−1、j≤−1相位、裁切端点和seam镜像预像；不冒充R41大型局部面积的独立第二实现。

`reflection_profiles_assembly_check.py`实际用Fraction组合15个operations JSON和轴误差JSON，exit_code=0，证据与范围见输出。单根式checker另为完全独立的有限bbox实现，与operations dyadic求轴核不同。

建议归档本人尾界理论/参数脚本/JSON、common144、单根式checker/JSON、本文、smoke checker/JSON、assembly checker/JSON。轴误差replay是作者JSON的重复副本，仅本地留证，不必重复入Git。R30/R37材料未改；没有提交或推送，正式main整合与远端核验由root处理。

复现接口：单根式checker支持`--catalog`（默认相邻common144）及可选`--seam`；小核smoke支持`--kernel-dir`（默认归档同computation）；组合checker要求`--profile-dir`和`--axis-error`。这些CLI不依赖作者clone的固定路径。


## 9. 冻结R46主稿最后绑定与双模板独审结论

最终完整审读主稿`proofs/high-ratio-reflection-edge-obstruction.md`，作者工作区文件13870 bytes，SHA256：

`9bf263dc2914990cd2643dd35f043aa29c5ac83f5ee817257f3b267c7ed9f8c8`。

主稿§§1–9通过。先前本review只独立组合R41；为覆盖冻结主稿的双模板结论，本次又独立读取R35全部14个新state、其已冻结bulk原记录（SHA256 24b1ddb373ba772588037f801a7b063e8fdddd170cee5e0defed5627169e7d74）及R35轴误差记录，准确复算全部区间。R35 bulk是明确复用而非本次重跑。最终共用独立checker为`computation/reflection_profiles_assembly_check.py`，替代归档中仅支持R41的早期本地checker；准确结果为`results/reflection-edge-kernel-independent/r035-assembly.json`与`r041-assembly.json`。

R35真轴leading上端约−9.105721697777488，F6上端约−53.293751299325；R41对应约−9.105721600526445、−53.2937507507865。两次实际Fraction运行都直接断言leading上端<−9、F6上端<−53，exit_code=0。故对两个固定模板每个整数N≥6，准确有

Jref<A/m−(9N−1)/(mN²)。

主稿中的稳定族上确界=A/m也正确：同一个准确面积多项式给Jref→A/m，而每个有限N≥6都在其下。这是族内上确界与不达到结论，不是一般C2或一般L¹不达。

## 10. 主稿§8有限L¹外扩上界的独立核对

设原子ν有n个正标签、总质量T，云团支持在各中心[-η,η]²。任一云团严格见证半径r捕获质量>4r²。取所有与该方形相交的云团标签I，它们全质量M_I至少为已捕获质量，故r<r_I=√(M_I/4)，各中心距接收中心≤r+η<r_I+η。因此云团严格halo包含于全部原子接收矩形各面外扩η后的并；即使原矩形为空、经外扩才激活也包含在此表示中。

原子子集r_I≤√T/2。若原矩形两边均正，外扩面积增量≤2η(W+H)+4η²≤4η√T+4η²。若原来一边非正、另一边可能为正，外扩后面积≤2η(√T+2η)，仍小于该统一界；两边非正则≤4η²。因此全部2^n−1子集逐盒付费确实给主稿上界，包括空盒/激活，不能只计原非空矩形。

η≤1时，新增halo面积≤4(2^n−1)η(√T+1)。按主稿η上界，此新增面积≤9N−1；除总质量T=2mN²，至多补回(9N−1)/(2mN²)，故严格保留原子亏损的一半。所选正有理η存在，也可用√T的有理向上界直接给有理η。这个上界极粗但方向、归一化与常数都正确；不排除任意尺度平滑后的其他族。

最终主稿与本review的范围一致：两个明确输入、规范N×N相位、一次真正最右反射、N≥6；N≤5及其他路线开放。本人冻结R30/R37七文件hash再次准确核对，全部未改。新独审材料未提交或推送；root统一整合。


## 11. 切相位补充的最终独审绑定

主稿冻结字节不改；另实际完整读取补充`proofs/high-ratio-reflection-phase-supplement.md`，SHA256：

`8730f0c3a6ce644b1ba0ad8262a8cde6cfcc3c88571ea884faf567f3fdc91e98`。

该补充明确采用本输入的充分条件C+4Px>30，避免只写周期长度而未展示切相位余量。本人另独立读取两冻结输入，以Fraction重新准确求得

xmax+4Px−30=424640887836625294630001/150000000000000000000000>0。

同时核对补充所列输入hash、canonical横相位与作者小checker/JSON的代码hash和准确分数；两输入均实际PASS。因为C>xmax，故对真轴及任何外移轴c≥C有

c_N−30=(N−5)Px+(c+4Px−30)>(N−5)Px。

左半接缝窗中的任何r<15严格见证，其原侧source满足x>c_N−30；canonical source cell编号≤N−6的点则x<(N−5)Px，故被排除。镜像侧取反射预像给同样限制，旧r<10目录亦被较宽R15 source-window涵盖。N≥5时下cut严格>0，因此确实看不到最左截断；与3Py>15、N≥6纵向分解合并后，主稿面积多项式、系数及全N≥6结论保持不变。

补充通过独审。没有重跑大型state，没有修改主稿、source kernel或任何R30/R37冻结文件；本次只追加本节及更新交付manifest的review条目。
