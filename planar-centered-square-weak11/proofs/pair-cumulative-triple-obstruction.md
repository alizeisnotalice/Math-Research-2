# 全部二阶Möbius累计项仍有严格三源strip欠账

日期：2026-11-15（Asia/Shanghai）。作者工作区 `pcw-sol-recursive-counterexample`。本稿只否定一个明确二阶累计模板：单源halo累计势加上全部pair的signed Möbius halo累计势，并仅在pair的源端点柱分摊该pair的总charge。它不否定任意两柱函数构成的势，也不声称任意pair储备机制都失败。第42轮文件冻结未改。

## 1. 模板和空中心的可核验需求

在阈值1、轴平行中心共同半径方形下，令A_I为源子集I的完整H_R在当前strip中的面积。单源累计势P₁满足
\[
 \Delta P_1=\sum_i A_i-\text{当前源柱charge}.
\]
对每个pair定义signed二阶密度
\[
 d_{ij}=A_{ij}-A_i-A_j.
\]
pair累计项满足\(\Delta P_{ij}=d_{ij}\)减去分配到当前柱的端点charge；每个pair的总charge可按任何固定方式分到两个源端点。即使pair总Möbius面积为负，也允许signed charge。本例当前柱没有任何源，所以这些端点charge全部为0。

对于N=3，模板在这一条边上的准确势差因此为
\[
 \Delta(P_1+P_{12}+P_{13}+P_{23})
 =A_{12}+A_{13}+A_{23}-A_1-A_2-A_3. \tag{1}
\]
中心真实质量也为0，故任意z的证书都必须令真实A₁₂₃不大于(1)。下面直接计算全部面积，不能仅从某点不被pair覆盖就推断(1)欠账。

## 2. 三源、完整尺度与strip

取h=1/2、T=8、R=2，三原子及各自纵向8周期lift为
\[
 p_1=(-1,0),\quad p_2=(1,3/2),\quad p_3=(2,0),\qquad m_i=4.
\]
当前strip为D₀=[0,1/2)×[0,8)。源柱索引是−2、2、4，当前柱0为空。k=ceil(R/h)+1=5；三个源在前后状态中均保留：状态长度是2k=10，前状态索引−5,…,4，右移后这些源位于−3、1、3。故此欠账不是端部源进入/离开导致。

不同纵周期的原子最小距离为8−3/2=13/2>2R，所以r<R不可能跨纵周期捕获。可在一个平移后的y周期内积分，也等于把这三个有限原子的halo沿全部y积分。单源半径1、pair半径√2、triple半径√3，全部小于R；全部7个非空源子集给完整halo。

如需全r>0，r∈[2,13/4)仍只能捕获同一周期源，总质量≤12，平均≤3/4；r≥13/4时逐纵轨道计数给
\[
 \mu(Q)/(4r^2)\le 3/(4r)+3/r^2\le87/169<1.
\]
因此本例H_R=H，无大半径遗漏。

## 3. 全部singleton和pair strip面积

三个singleton接收方块边长2。在当前strip内，仅p₂的方块有正面积：
\[
 A_1=0,\quad A_2=1,\quad A_3=0. \tag{2}
\]
令a=√2−1。pair(1,2)的共同接收矩形为
\[
 E_{12}=(1-\sqrt2,\sqrt2-1)\times(3/2-\sqrt2,\sqrt2).
\]
它在p₂ singleton外留下
\(x\in(0,a),\ y\in(3/2-\sqrt2,1/2)\)，面积a²=3−2√2。pair(1,3)水平跨度3>2√2，矩形为空；pair(2,3)共同矩形左端2−√2>1/2，不进入当前strip。因此
\[
 A_{12}=4-2\sqrt2,\quad A_{13}=0,\quad A_{23}=1. \tag{3}
\]
所以(1)的二阶预测恰为
\[
 A^{(2)}=4-2\sqrt2. \tag{4}
\]
当前strip中的singleton/pair halo并也恰为p₂ singleton并E₁₂，面积同为(4)。因此本例没有被pair Möbius负项或pair之间重叠抵消掩盖：purepair union包络和signed二阶累计模板在此给出同一个准确预测，尽管它们一般不是同一对象。

## 4. 完整triple新增与准确残差

三源共同矩形为
\[
 E_{123}=(2-\sqrt3,\sqrt3-1)\times(3/2-\sqrt3,\sqrt3).
\]
在strip中，去掉p₂ singleton后留下一个面积
\[
 Q=(\sqrt3-3/2)(\sqrt3-1)
\]
的矩形。该矩形与前述pair新增区域的交面积为
\[
 O=(\sqrt2+\sqrt3-3)(\sqrt2-1).
\]
全部相对端点次序由√2∈(7/5,3/2)、√3∈(17/10,7/4)及相应平方比较确定；有限条接收边界线不影响面积。故完整strip面积及二阶欠账准确为
\[
 A_{123}=\frac72+2\sqrt2-\frac{3\sqrt3}{2}-\sqrt6,
\]
\[
 \boxed{A_{123}-A^{(2)}
 =Q-O=-\frac12+4\sqrt2-\frac{3\sqrt3}{2}-\sqrt6
 >\frac{211}{1980}>\frac1{10}.} \tag{5}
\]
严格有理下界使用√2>140/99、√3<26/15、√6<49/20；三个比较分别平方即可证明，代入得到560/99−111/20=211/1980。

由于当前中心柱质量为0，任何z均无法补齐(5)。在本例这同时是完整halo超出全部single/pair halo并的真实面积，因此purepair union版本同样欠账；结论依据完整并面积，而非仅有点态纯三源见证。

## 5. 显式严格共同激活矩形

作为独立几何检查，取
\[
 U=(7/16,1/2)\times(-1/32,1/32),\quad |U|=1/256.
\]
固定r=13/8捕获三源，其平均为192/169>1。对任何x∈U，p₁、p₂、p₃各自捕获距离分别大于23/16、47/32、3/2，均大于单源临界半径1。pair(1,2)捕获p₁所需半径>23/16>√2；pair(2,3)捕获p₃所需半径>3/2>√2；pair(1,3)原水平跨度3也需r≥3/2>√2。因此U在全部singleton/pair halo之外，且属于真正共同半径三源halo。

这个正面积见证单独只否定single/pair集合覆盖；式(2)–(5)另证明没有面积抵消，才升级为所限定累计势模板的严格失败。

## 6. 状态要求与证据范围

单源累计势能按√m改变有效尺度。要准确容纳某个pair，其水平捕获范围至少须考虑
\(\min\{2R,\sqrt{m_i+m_j}\}\)：若两源水平跨度不小于该量，则其pair共同矩形为空；反向还须检查纵向跨度。这个质量适应范围能逃过固定物理距离的截断，但本例说明即使全部pair都已精确累计，仍有真实三阶缺额。需要新的三源储备、允许非端点收费的重分配，或另一种整体几何证书；这里未证明哪条机制足以处理一般输入。

脚本 `computation/check_pair_cumulative_triple_obstruction.py` 实际运行Python 3.15.0、标准库Fraction和整数isqrt，单进程，exit_code=0，无随机种子、无浮点容差。对全部7个非空源配置再完整枚举其全部非空源子集，共19个原始子集；用D=2^44的定向有理平方根内外包围及独立纵条带并面积算法，严格确认(5)与超出purepair union的面积均>1/10；并检查上述严格有理见证及尾界。没有求解大图，也没有把有限例当作一般pair模板不可能定理。

运行命令：

```sh
python3 -I -B planar-centered-square-weak11/computation/check_pair_cumulative_triple_obstruction.py
```

脚本SHA-256：`f5b3e427d75868790ec86820b7212cc9254cbc9d891240108f08cd257f54654e`。结果 `computation/results/pair-cumulative-triple-obstruction/exact-check.json` SHA-256：`7be384ff418740a06bf08bdc1038da243ea8b8d5a6798a18aecf310f27e8870f`。本稿未提交或推送。
