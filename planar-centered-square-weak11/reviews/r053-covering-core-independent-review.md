# R53覆盖核心独立审查：九格全质量证书与定量bbox收缩

2026-11-15。最终判定：下述明确范围内通过。九格全部闭质量单纯形已经严格闭合，得到 B<33/10；任意大bbox输入的一般C₂仍开放。作者生成器、作者重放器、原始证书均未由本审查修改；我的两份检查程序在独立clone，实际运行由operations协调。

## 一、闭单纯形覆盖与根式外矩形

审查对象为全部 w_i≥0、Σw_i=1，包括所有零质量面。原始177节点证书与加强32251节点证书采用同一严格机制，检查程序从根盒[0,1]^9重建。

每个节点的下界 max(l_i,1−Σ_(j≠i)h_j) 与上界 min(h_i,1−Σ_(j≠i)l_j) 都由更新前的父盒计算；对任何满足总和1的可行向量成立，因此收紧不删除可行点。每次分支用两个闭半盒，分割cut由两边共享，故不会漏掉单纯形边界、分支面或零质量。访问全部节点恰一次，根到叶深度逐项相符；叶排列是0,…,8的完整置换。

叶前缀上界 u_k=min(前缀上端和,1−后缀下端和) 对整个该叶成立。检查器重新计算有理u_k，然后只用整数平方不等式检查矩形所用半径H_k/D是最小向上根界；不调用生成器的isqrt/ceil例程。坐标格端点、矩形四边、D=2^40全部独立重建。保留全部九格，即使w_i=0；这扩大目标并保持向上界方向。实际跳过零质量格后，后续非空格累计质量不变。

面积核验采用**所有511个非空矩形子集的整数容斥公式**。每个交集面积由 max下端、min上端及正部直接算出；不调用作者的横扫或端点网格helper。逐叶容斥整数与证书面积一致。取得的统一加强上界为

    B = 1994727599733432151936139 / 604462909807314587353088
      ≈ 3.299999995647862 < 33/10.

B是证书上界，不是自由顺序最优常数或全质量最坏值。与7/2的准确余量为

    120892584592168903799669 / 604462909807314587353088.

## 二、真实proxy与真实halo映射

在边长√M的ambient内，用quarter/middle/quarter坐标带给源质量唯一归格；闭格只用于外容器，cut处不重复记质量。按证书排列整格插入有限原子时，格内每个proxy有效质量≤该格完成前缀质量Ms_k、中心属于相应闭格。因此真实proxy并包含于√M缩放的证书容器，得到≤BM。格内原子数、单原子质量、位置及eligibility rank都不受限制。

还有直接适用于任意有限正Borel测度（含L¹密度）的同一局部定理，无需先把Borel测度近似成原子。给定真实严格见证Q(x,r)，取证书顺序中最后一个被Q捕获正质量的格k。则捕获总质量≤Ms_k，故r<√(Ms_k)/2；该格与Q有非空交点，得到x在该格闭bbox加该严格半径的容器内。对全部见证取并，H(μ)包含于同一个证书容器。因此只要支撑可装入边√M方块，即有|H(μ)|≤BM。门槛λ>0时把μ除以λ，条件为支撑可装入边√(M/λ)方块，结论为≤BM/λ。

## 三、定量bbox收缩证明独立核验

结构作者给固定不交源标签块C_i、质量m_i、真实bbox B_i、代表点c_i∈B_i，以及最大位移δ_i。D_J=∩_(i∈J)(B_i+开方形半径√m_J/2)。真实严格见证实际遇到的块集合J使r<√m_J/2，故其中心属于D_J；完整块质量用于保守外包，并未假设见证捕获整个块。

代表点完整子集公式给H(μ_c)=∪_J E_J(√m_J/2)，且E_J⊂D_J⊂E_J(√m_J/2+δ_J)。因此|D|−|H(μ_c)|≥0，并由各E的几何外差之和支配。这里求和的是外差面积，不是将重复m_J当成新源质量。若保留全部交叠，直接算|D|−|H(μ_c)|更紧。这只传递真实H，不承诺原插入proxyU在压缩后不变。

饱和K(J)的四个端点阈值条件均保持原极值，因为J⊂K且原极值有J中的达到点；质量只增。故D_J⊂D_K，且K幂等。至多min(2^n−1,n^4)个不同状态的断言成立；不能擅加bbox有共同交点。δ可以超过eligibility余量，代表点可能使eligible子集翻转，证明仍不需要任何最小质量或余量。

误差公式为两条正部边长的乘积之差。若原两边非负，展开给e_J≤4δ_J√m_J+4δ_J²；若某边为负，扩张后该边≤2δ_J，另一边≤√m_J+2δ_J，得到更小的界。零边和原空集均已包括。定理可以接受任意高rank、大bbox的固定分割，但其误差不自动小；n随细分增长时也不能由δ→0擅称总误差趋零。

## 四、已实际核验的适用范围扩展

新的局部定理有简单的质量放大推论。令d为整个支撑bbox最大边长，a=max(1,d²/M)。由μ≤aμ及|supp μ|的几何不变，放大测度aμ质量为aM≥d²，故

    |H(μ)| ≤ B max(M,d²).

因此d²≤(35/33)M时严格有|H(μ)|<3.5M。有限原子输入若每个实际eligible超图组件分别满足此条件，则组件halo分解与面积次可加性给同样的一般输入条件类上界；组件不必整体eligible。对完全不满足几何条件的组件，此推论没有统一常数。

另一保守但严格的扩展：若至多9个源块可压到边√M方块内的代表点，且各位移δ_i≤√M/20000，则q≤511，收缩误差≤4·511·(1/20000+1/20000²)M<B与3.5之间的余量。于是原输入也有<3.5M。这个数字只是显式充分条件，不是最佳位移容忍度；一般大bbox组件未必能分成9个这样的紧块。这个紧块例的原始bbox也满足d_original≤(1+2/20000)√M，因此已被上一段质量放大推论更强地覆盖；它展示收缩误差的显式支付，不声称增加一个超出质量放大条件的新类。更实际的递推应使用上述bbox叶和准确接收遮罩，不能把有限单纯形定理直接推广成C₂≤3.3或3.5。

## 五、实际运行与版本

原始独立检查（同为Python3.13.18）：177节点、89叶、801个向外半径、45479个容斥子集。加强检查：32251节点、16126叶、145134个向外半径、8240386个容斥子集，实际11.5614秒，Python3.13.18，单核，退出0；标准库整数/Fraction，无随机数或浮点判断。另有576个有理标量分支测试支持收缩误差公式，有限测试不代替上一节的解析证明。加强wrapper还准确核对35/33的bbox质量放大及1/20000位移预算。

代码为 computation/review_r053_mass_simplex.py 及 computation/review_r053_mass_simplex_33over10.py；第一份在原输出后字节冻结，第二份导入其run函数，额外处理gzip与33/10目标。输出分别为 computation/results/round053/independent-mass-simplex-review.json 与 independent-mass-simplex-33over10-review.json，后者SHA为2471eb2957d16877d5b56706f1e7ad355145ee2e61c6d03f4c832d11d3f786ee。运行在operations工作区，现复制到本clone的输出字节完全相同。

复现可指定CodexTemp输出，不覆盖归档记录：

    python3 -B computation/review_r053_mass_simplex_33over10.py --certificate computation/results/r053-mass-adversary/all-mass-box-certificate-33over10.json.gz --out /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/r053-independent-replay.json

本审查没有重跑作者Dirichlet/DE探索，没有认证其全局最优，也没有重跑旧一般核上界。主目标仍需要任意输入的无重复质量几何分解或接收预算；以上成果是完整闭合的局部定理及明确损失接口。

## 最终来源绑定

- `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-cover-counterexample/planar-centered-square-weak11/proofs/r053-nine-cell-mass-simplex-certificate-33over10.md`：`f7e9223403939788d940ebc1cfd47d69ad772eb44da73e42cdeb08d6f80fd784`。
- `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-cover-counterexample/planar-centered-square-weak11/computation/results/r053-mass-adversary/all-mass-box-certificate-33over10.json.gz`：`e6dd2d93ef753857e6645a6eeb12a4e23755c762f21d8a7bf7722c377ef9e2a2`。
- `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-cover-structure/planar-centered-square-weak11/proofs/r053-quantitative-bbox-contraction.md`：`c7aafe8a211a39e15b3f5980f6bac0b03f4c444bc5cfef2a0a2a5aee8608ac10`。
