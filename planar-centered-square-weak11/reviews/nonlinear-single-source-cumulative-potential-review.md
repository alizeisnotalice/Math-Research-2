# 第42轮独审：非线性单源累计势与hardcore pair证书

审者：sol_recursive_geometry；日期：2026-11-15。结论：**§§1–5数学通过独立解析与小型精确计算审查，无必须修复项**。不把hardcore受限类证书解释为一般有限图层级上界。

源文件：`/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-recursive-counterexample/planar-centered-square-weak11/proofs/nonlinear-single-source-cumulative-potential.md`。

最终绑定SHA-256：`c2c2cab8af8a3b1964b3fafe15f95760c6a71b5125990714384930c4d386d6eb`。基础版本 `9796de7d6c69bc257aca0a6e88f3d1c267c11a2cc84e4d84bb870d53dfd943d0` 的§§1–5已完整审读；最终回读新增(6a)(6b)单柱弱周期尾界及实际范围和计算段后，绑定以上最终版本。

## 1. 累计φ的质量依赖、差分及边界

在T>2R下，r<R不可能同时捕获同一柱两份纵向lift。孤立原子的完整capped halo确为半宽t=min(R,√m/2)的开方块，纵向一周期的两半总高2t。a_j公式正确，按平移条带分割求和得A=4t²≤m；m=0无歧义。

φ_j=Σ_{l>j}a_l−A1_{j<0}精确给φ_{j−1}−φ_j=a_j−A1_{j=0}。顶点slot j在下一状态读取物理柱j+1，因此P₁差分中每物理柱的系数是φ_{j−1}−φ_j，符号正确。

边界：k=ceil(R/h)+1使a_j在j≤−k、j≥k都为0；因φ的累计尾和，此时φ_{−k−1}=φ_{−k}=φ_{k−1}=φ_k=0。故进入/离开状态时未遗漏源或成本。作者正文所写两处边界零值足够，并可如上补全直接用于差分式。

a_j=a_{1−j}，故j≥1部分总面积为A/2，j≤0部分亦A/2。对j≥0，φ是非负尾和≤A/2；对j<0，φ=−Σ_{l≤j}a_l，绝对值≤A/2。由此|φ_j|≤A/2≤m/2。若一个窗至多一根非零柱，a=S，差分=S−A(m₀)，证书a−m₀≤ΔP₁只有slack m₀−A(m₀)，未增加虚构质量。

m_l=4l²h²时a_j=2lh²1_{1−l≤j≤l}。累计式的负半轴端点是(l+j)_+而非(l+j+1)_+，作者公式正确。z=1同时支付全部这些single-column测试，所需非线性在m上，未与线性系数障碍矛盾。

## 2. 弱周期条件与实际第32轮字母表（新增段落）

不必T>2R：对一柱每周期仅一个原子，若m≤T²/2，则r≥T/2的方形至多捕获2r/T+1份lift，故平均≤m/(2Tr)+m/(4r²)≤2m/T²≤1。最后允许等号，因为水平集严格>1。

r<T/2最多一份lift，又√m/2≤T/(2√2)<T/2。因此完整halo仍由各孤立singleton方块构成，capped版本a_j、A和φ公式完全保留。

对第32轮n≥1，m∈W_n有m≤2^n≤T_n²/2，且R_n>√m/2；所以A(m)=m，P₁确实能直接支付实际第32轮全质量字母表的**单原子柱**single-column窗。仍不包含同一柱多点profile。

随后(6)的m_l闭式在§1原域T>2R或(6a)替代域m_l≤T²/2内调用；不能把单独lh<R、2lh<T当作排除多lift的充分条件。作者保留原域并明确给出替代域，本审按此前提核验；实际第32轮(6b)满足替代域。

## 3. 空中心真实mixed缺额

h=1，源在−1和1，各m=9/4。singleton半宽3/4：左源不交中心strip，右源交宽3/4、高3/2，S=9/8。

pair门槛t=3√2/4>1；共同接收矩形E=(1−t,t−1)×(−t,t)。其水平支撑在两singleton的空隙内，因为t−1<1/4，故与singleton内部不交。中心strip保留该pair矩形的右半，准确B=2t(t−1)=9/4−3√2/2。√2<17/12使B>1/8，严格余量正确。

全isolated pair的完整面积为2m+2B=9−3√2，中心strip完整面积9/8+B=27/8−3√2/2；三个非空源子集已全部覆盖。m₀=0且ΔP₁=S，因此任意z都会留下正缺额B。单纯提高z无法修复空中心边。

有理开矩形U与r=33/32见证也正确：两原子都被同一个方形捕获，平均128/121>1；每个singleton的捕获距离≥31/32>3/4，所以U严格在两独立halo外。没有混用不同接收点或不同轴半径。

## 4. 任意h的hardcore完整几何

约束c_jc_{j+1}=0。r<2h时三源横向跨度≥4h，无法共同捕获；不同纵向lift相隔8h，更不能共捕获。pair间距≥3h的最小半径≥3h/2>3√2h/4，故不激活。只有间距恰2h的pair有效。

singleton之间间距≥2h、半宽3h/4，内部不交。有效pair矩形半横宽t−h；它到任一源中心的横距至少h，而singleton到其源仅延伸3h/4；t<5h/4，故pair与全部singleton不交。不同有效pair中心间距至少2h，各宽2(t−h)<h/2，亦不交。纵向周期lift不重叠。

所以完整halo严格等于这些互不交的singleton与有效pair矩形并。中心strip只收到源0、1的各半singleton及pair(−1,1)、(0,2)的各半矩形，恰得a=m(c₀+c₁)/2+B(c₋₁c₁+c₀c₂)，没有漏计相邻mixed或多源见证。

## 5. 势差、零中心、周期接缝与全尺度tail

当前质量下φ仅j=0非零且为m/2。对p=m c₀/2+B c₋₂c₀，右移后pair项为B c₋₁c₁，故

    Δp−(a−z* m c₀)=B c₀(2−c₋₂−c₂)≥0.

特别c₀=0时slack准确0，真实pair(−1,1)收益由新状态pair项准确支付，没有凭空给空柱记质量。

长度6无相邻1字共21个，长度7共34个；有向重叠把所有周期接缝的硬核约束纳入。周期1010逐边等号；每源左右距离2h的pair均存在，所以其比值等于z*=1+2B/m=3−4√2/3。沿任意闭游走求势差得上界，1010给匹配下界；结论对所有h>0按h²缩放成立，只限该binary/hardcore/single-y0类。

全尺度tail独立成立：r≥2h时每行横向源间距≥2h，源数≤r/h+1；每列lift数≤r/(4h)+1，乘以m并除4r²得9/64+45/(64u)+9/(16u²)，u=r/h≥2。它递减，u=2时81/128<1，连r=2h也严格排除，故开cap H_R与全尺度H确相等。未把r<2h几何直接当作大半径结论。

## 6. 一般profile条件接口

用孤立柱实际完整capped strip面积a_j(c)、A(c)=Σa_j(c)替代单原子公式，累计差分仍成立，边界只依赖横向最大半径R。此时A(c)≤m(c)不再自动成立；作者明确要求κ_col及z≥κ_col，纵向双原子反例可使κ_col>1。

D(w)=a−Σa_j可正可负，加入mixed势后的必要充分逐边预算D≤ΔP_mix+z m₀−A(c₀)只是接口，不是已证一般上界。作者最后范围表述正确。

## 7. 独立实际计算与证据范围

新写并实际运行 `computation/check_nonlinear_cumulative_pair_review.py`；Python3.15.0，标准库Fraction及独立Q(√2)比较/运算，单进程exit_code=0，无随机种子、无浮点容差。

没有调用作者PASS或复用其几何输出：本脚本从每个源子集的外接区间与singleton/pair门槛生成接收矩形，在lift y=0、T上裁切，再以独立x扫描+y区间并求完整面积。≥3源和跨lift由上述解析cap条件排除。

实际核对：21顶点/34合法窗；h=1、5/7两个尺度共68次完整几何/势差，含50次空中心边；7种质量根含3个cap激活/边界例，91条累计恒等式；48条弱周期尾界与第32轮n=1,...,9最大字母质量条件；两个源的中心strip与全平面准确面积；tail81/128与有理严格见证128/121。

结果：`computation/results/nonlinear_cumulative_pair_review_exact.json`。
复现：`python3 planar-centered-square-weak11/computation/check_nonlinear_cumulative_pair_review.py`。

本审没有重新跑作者4192单源测试/102次几何，将其视为作者证据，未冒充本审执行。有限计算辅助解析独审，不代替任意h、任意质量与全部窗的一般证明。

未修改作者文件或第32轮冻结文件，未提交、未推送。一般多profile、多质量及多源mixed统一预算仍开放。
