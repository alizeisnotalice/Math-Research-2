# 八原子全权重有限集合质量下界：独立256bit复验

日期：2026-11-15。第51轮独立附审；独立工作区`pcw-sol-kernel-certificates`。结论：作者187矩形有限测试集合证书通过独立积分和几何复验，严格给出所有θ∈[0,1]的I_θ>40/19。与已证明的R/(2m_θ)>19/10结合，整个同L1半径轴/对角八原子混合族的合法标量majorant成本严格>4。本结论不是C₂下界。一般正部势的生成元分布合法性由主组另一份独立证明承担，本附审不把该接口当作本次重新认证的内容；必要成本排除本身不需要核合法性。

## 1. 真实来源和输入

只读作者工作区`pcw-astra-kernel/planar-centered-square-weak11`：

* `computation/cauchy_eight_atom_mass_floor.py`，SHA-256 `eb07a10a8a933cbdf155226f8125be6ffee3457df04c43abb484e36e8e3932a5`。
* `computation/results/cauchy-eight-atom-mass-floor.json`，SHA-256 `a8aa37d3c78e164f8e8c0510e06e5d474093723a49e0f0b80d45b08c8258286c`。

本次没有取原JSON的区间值参与任何算术，也没有导入原signed H程序。仅提取187个矩形的准确有理端点及4/8权重，另存`computation/inputs/cauchy-eight-atom-mass-floor/rectangles.json`，使复现不依赖作者的数值输出。N=8和中点θ=1/4只说明作者如何找到S；完整证明允许任意可积S，不需要重跑选择过程或假定S内逐点正性。

规范固定R=1、U=1/(|x|+|y|)，axis四点(±1,0),(0,±1)，diag四点(±1/2,±1/2)，每端测度总质量1。I_θ=∫[U-U*(θν_axis+(1-θ)ν_diag)]_+。

## 2. 独立积分实现

独立程序`computation/cauchy_eight_atom_mass_floor_independent.py`对每个矩形和每个源，先在x=源横坐标、y=源纵坐标处分割。每小块的两绝对距离各单调，将它们实际映成非负距离盒[a,b]×[c,d]。再用

    F(t)=tlogt-t，F(0)=0，F''(t)=1/t，
    ∫∫1/(u+v)=F(b+d)-F(a+d)-F(b+c)+F(a+c)。

该式在源落角点时仍由可积极限成立。它在实现结构上不同于作者的全平面带符号H；没有借用其corner evaluator。所有端点是准确有理数；Arb只承担log和算术的向外误差。逐矩形求原点势积分，再分别扣去4个源积分的平均，最后乘真实轨道权重并累加。

## 3. 集合几何及权重

187个原矩形分成25个中央第一象限单位方形内矩形和162个x>1、0≤y≤1的右臂矩形。前者分别对两轴反射4倍；后者先取坐标交换或不交换，再分别反射，展开8倍。

本次实际生成1396个展开矩形。所有有理端点乘8后精确成为整数；逐对准确检查973710对，没有一对内部具有正面积重叠。展开面积为2709/16，也与原权重面积和准确相等。因此4和8不是抽象乘数，而是一个真实有限集合S的无重复面积表示；共享边角为零测。

W_axis和W_diag分别在反射和坐标交换下不变，所以每个展开像上的积分等于原矩形积分。中央区与各臂在内部不交，其不同反射扇区也不交。实际逐对检查进一步排除了选中矩形之间的漏计重叠。

本证书是对有限S的下界计算，无需覆盖整个正区或给无限尾上界。S外可能有无限正尾，舍去它只减小正部总积分；没有把x≤32宣称为核支撑。

## 4. 实际256bit结果

实际解释器Python3.12.13，python-flint0.9.0，ctx.prec=256，单进程，exit code 0，记录耗时约0.120秒。独立得到

    V_axis=∫_S W_axis
      ∈[2.1562814549962867652703622591711772800037827192682319364389556533405719
             +/-5.62e-71]，
    V_diag=∫_S W_diag
      ∈[2.11590746225515448882397803939726837532992241406826722140029952339546001
             +/-9.35e-72]。

这里的十进制仅便于阅读。程序把真正Arb向外端点转成准确dyadic有理数，并用Fraction直接认证两下端严格>40/19；每个矩形的两端积分和最终端点均保存在实际输出中。

对每个实θ∈[0,1]，点态正部不等式给

    I_θ≥∫_S W_θ=θV_axis+(1-θ)V_diag>40/19。

没有把I_θ的凸性误用成端点下界插值。这里线性的是固定集合上的势差积分，严格下界对整个连续θ区间成立。

## 5. 成本结论及可复现范围

R51作者主稿`pcw-astra-kernel/planar-centered-square-weak11/proofs/cauchy-eight-atom-positive-potential-floor.md`已完整只读审查，绑定SHA-256 `b2651f6fd952623977c8853e3018262ea6dd28539f66ed81ebe4e75a3a1cde60`。其证明仅使用轴顶点的必要值v_θ(R)>0及a≥1/v_θ，不把顶点值当作整D最小值。对R≥2，v_θ≤m_axis，P(R)>0给R/(2v_θ)>19/10；1<R<2，两端顶点值均≤R/4，故该因子≥2。公式、差值、严格方向及P的正平方分解均通过核对。

因此R51排除结论不依赖较强的整D准确最小值定理。该定理另保留在`proofs/cauchy-eight-atom-diamond-minimum.md`供可行性与优化使用，导数和分区证明亦已独立完成。用必要顶点值和本次质量证书已经得到

    J_θ> (19/10)(40/19)=4，

覆盖所有R>1、θ∈[0,1]的该同半径族及其两个四原子端点。不同轴/对角半径、其他原子位置、非正权重和一般势差核不在本结论范围。

复现从课题目录运行：

    /Volumes/data/CodexProjects/T/planar-centered-square-weak11/computation/kernel-baseline/venv/bin/python -B \
      computation/cauchy_eight_atom_mass_floor_independent.py

实际结果：`computation/results/cauchy-eight-atom-mass-floor-independent/checks-256bit.json`。解释器是用户已有环境；没有原地运行旧T脚本。未安装软件、未运行MATLAB/LP/大积分、未改作者或main/R33/R38冻结材料；未提交推送，由root白名单归档。一般C₂精确值仍开放。
