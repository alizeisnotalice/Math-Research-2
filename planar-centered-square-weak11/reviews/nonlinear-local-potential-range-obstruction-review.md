# 第50轮独审：非线性单柱势的物理作用范围必要增长

审者：sol_recursive_geometry；2026-11-15。结论：**全稿主命题及新增§5式(8)通过，无必须修复项**。

冻结源：`pcw-astra-recursive/planar-centered-square-weak11/proofs/nonlinear-local-potential-range-obstruction.md`，SHA-256 `e0038121c05263103f64b421ee13f92c6ebaf7d9a4c00a64cb4338eb6319c781`。

本审先读旧版 `4d7a1671ac56ee8330f22ce007a79ade7ba814afcbbf4e06f900b916f49990e9`，后重新回读§5质量量词修订并绑定以上最终稿。最大单点质量2^n不能混称一般柱profile最大总质量；后者有2^{2n}个格点，总质量可达2^{3n}。最终稿已明确前者，证明全程只用单原子profile，适用范围准确。

## 1. 对独立一般证明的复核

本审此前已独立证明§§1–4同一命题，附件 `proofs/fixed-physical-range-potential-obstruction.md` SHA `224227a560a7199a06f6d284f83c2d1cd83810667c2b3b0352d4d7d0de2aaa2d`。

回读主稿确认：φ_j(0)=0已规范化；固定同一质量后α_j=φ_j(e_m)/m只是函数取值，没有把非线性势全局线性化。稀疏窗两端I均为0，单源迫β≥0，Σβ=z。每dyadic环带全部q余数类必须都求和，不能只用单个类而漏掉q；主稿精确完成该步。每类K-subset出现概率K/M是组合恒等式，包含K=M端点。各环带半开故无重复。

共同临界半径2r_s、纵向两份lift总高4r_s及真实cap内严格见证都正确；只用内界，不需要tail。测试质量m≤2^n且为整数，W_n量化无损；窗口与周期条件均由同一质量约束保证。I的稀疏消失前提覆盖任意固定物理span的pure相互作用，不限pair。

本审此前真实运行 `computation/check_fixed_physical_range_potential_obstruction.py`，Python3.15.0 Fraction，单进程exit0：474个实际K-subset几何、36次组合平均、45个all-n质量选择、7个slot常数规范化。结果 `computation/results/fixed_physical_range_potential_obstruction_exact.json`。本次明确复用，不将历史运行冒称本次重跑。

## 2. 新增范围必要条件的完整代数

对第n层范围D_n，取整数L=floor(D_n)+1，满足L>D_n、L≤D_n+1。

若16L²>2^n，则L>2^{n/2}/4≥2^{n/2−4z}/4（用z≥0），故式(8)直接成立。

否则，令p≥0是满足16L²4^p≤2^n的最大整数。主命题给p+1≤4z；最大性给严格16L²4^{p+1}>2^n。由于4^{p+1}≤4^{4z}=2^{8z}，得到

    16L²2^{8z}>2^n,
    L>2^{n/2−4z}/4,
    D_n+1>2^{n/2−4z}/4.

严格号来自p的最大性而非连续近似，指数−4z及外部1/4均正确。若z小到无法满足p+1≤4z，则该证书前提本身不存在，必要条件当然仍成立。

本轮另在新脚本 `computation/check_pair_cumulative_triple_review.py` 中对4种L、n=1,...,24共96个分支/平方形式作精确有理算术检查；这是代数辅助，不替代上述任意实z证明。

对固定成本z，该条件确实迫D_n与2^{n/2}同数量级下界；不证明按此范围即可构造证书，不排除范围随实际捕获尺度或质量增长，也不排除每层完整势。最大一般柱质量的更大上限未被本测试使用。

未修改作者文件或第32/42轮冻结材料，未提交、未推送。
