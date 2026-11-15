# 第42轮独审：线性柱质量加纯相互作用势的障碍

审者：sol_recursive_geometry；2026-11-15。结论：**所审全部命题通过，无必须修复项**。

冻结源：`/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-recursive/planar-centered-square-weak11/proofs/linear-mass-interaction-potential-obstruction.md`。
SHA-256：`3faa10f5c216493d993da0c28c088131de9e8495271039380222a5bd5d754c39`。

## 核验

1. 顶点slots按−k,...,k−1编号，右移后同slot读取物理柱j+1；故β_j=α_{j−1}−α_j+z1_{j=0}，端点α_{−k−1}=α_k=0，Σβ=z。相互作用仅在单柱测试两端消失，作者明确未对一般多柱转移删ΔI，符号及量词正确。
2. 测试每周期唯一y0原子，质量4l²h²，临界半径lh。列l与1−l各在中心strip产生至少2lh²面积；y0、yT两份lift的各半段都计入。T>2lh仅保证这些singleton矩形不重叠，不需要排除一切更大半径mixed：任何额外面积都增强必要条件。严格接收边界的零面积线未冒充正面积见证。
3. 两组测试位置不交；每个其余位置用任意合法正质量单柱测试得到β_j≥0，故z≥H_L。每项是独立输入，而不是一个多柱窗，纯相互作用无机会支付。
4. 第32轮层n≥2取L=2^{n−2}，m_l≤1/4，且m_l/step=4l²2^{2n}整数，所有cap/周期/位置条件均合法。按dyadic blocks求H_{2^{n−2}}≥n/2成立。该发散仅为受限势函数模板的成本，不是一般有限图最优值发散。
5. §5限定每非空柱每周期只有单y0原子，范围明确。R=1<T/2=2使不同lift不能共捕获；t_l=l/64≤1/2，准确a_j/m_l=1/(2l)当1−l≤j≤l，其余0。因此对每位置取最小可覆盖l即得文中β包络，和为H32。递归恢复α的末端条件正由z=Σβ保证。这里仅认证全部一柱测试，作者未宣称覆盖一般多柱边。H32=586061125622639/144403552893600>4（分子减4倍分母为8446914048239>0）。

## 实际复用证据与限制

本审此前已独立建立同一一般命题，附件 `proofs/singleton-column-potential-obstruction.md` SHA `aa879d256235cb73b21ab9b819143a72b7ba6a96a3494b8b63efab8d128c04cc`。
实际已运行 `computation/check_singleton_column_potential_obstruction.py`，Python3.15.0 Fraction、单进程exit0，n=2,...,9共510个两侧single-column测试及β望远镜边界符号。脚本SHA `1c40c054705f745db9397a5612f34762931eefddfeb38488e6e5a85c96d1dc9b`；结果 `computation/results/singleton_column_potential_obstruction_exact.json` SHA `72fc8db88f10d91c626d35c1f3cf8c2da3595da2c36f572f32f2619c10c70847`。

此次明确复用上述真实运行，未冒称重新执行作者4192项；§5包络由本审解析比较独立确认。源版本末尾的作者4192项运行仅属作者证据，未在本审重新跑其脚本。

非线性single-column项或能在单柱状态取非零值的几何势均未排除。未修改作者文件或第32轮冻结材料，未提交、未推送。
