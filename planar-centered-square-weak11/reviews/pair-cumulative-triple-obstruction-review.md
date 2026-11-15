# 第50轮独审：端点收费二阶累计的完整三源strip欠账

审者：sol_recursive_geometry；2026-11-15。结论：**完整几何、全部strip面积、限定模板失败及全尺度tail通过，无必须修复项**。

冻结源：`pcw-sol-recursive-counterexample/planar-centered-square-weak11/proofs/pair-cumulative-triple-obstruction.md`，SHA-256 `79c0768da77029d38af3254cfe72466bb27c228a788fe9e1b953da4fea8490c2`。

本审未依赖作者PASS或复用作者矩形并代码/输出。以下独立解析复核和精确代数计算只针对这一个完整三源配置及明确累计模板，不推广为任意两柱势不可能。

## 1. 量词、图窗和周期lift

h=1/2、T=8、R=2，源(-1,0)、(1,3/2)、(2,0)，每份质量4；strip=[0,1/2)×[0,8)。所有方形轴平行、中心、共同半径，阈值严格>1，capped半径0<r<2。

源柱索引−2、2、4，k=5，前状态slots−5,...,4都包含它们；后状态右移一格索引−3、1、3仍都包含。当前柱无源，真实中心成本为0，亦没有源进入/离开而造成的伪缺额。

不同周期标签原子最小纵距13/2>2R=4，故capped方形不能混合两个周期标签。每个同周期子集的临界半径为sqrt(|I|)，均<R。所有7个非空I及其所有非空子集，共19个原始目录项，确实覆盖完整halo。纵向有限halo平移8后不会互相覆盖；其[0,8)面积等于一个有限三源halo在全y的面积。

特别：triple接收矩形有一块y<0，它在[0,8)内由y=8的整组lift贡献。不能忽略这块，否则总strip残差会算错。独立脚本明确按lift−8、0、8裁切，并确认只有必要lift有非空片段。

## 2. 全部singleton与pair面积

单源临界半径1。源1方块右端恰x=0、源3方块左端x=1，都无strip正面积；源2交宽1/2、高2，所以(A1,A2,A3)=(0,1,0)。临界边线不计面积。

pair门槛sqrt2。E12的strip部分横宽a=sqrt2−1，纵范围(3/2−sqrt2,sqrt2)。在源2 singleton外只新增底部一块，纵高同为a，因此A12=1+a²=4−2sqrt2。

pair13横跨度3>2sqrt2，矩形为空；其两个singleton也都不进strip，所以A13=0。E23左端2−sqrt2>1/2，不进strip；源2仍贡献1，所以A23=1。

因此signed二阶预测Σpairs Aij−Σsingles Ai=4−2sqrt2。所有single/pair halo的集合并也恰是源2 singleton并E12，面积同为4−2sqrt2。两个对象在本例相等由完整目录确认，而不是一般恒等式。

## 3. triple矩形、重叠与整个strip残差

共同门槛sqrt3；E123横向(2−sqrt3,sqrt3−1)、纵向(3/2−sqrt3,sqrt3)。裁入strip并去掉源2 singleton后，剩余矩形宽sqrt3−3/2、高sqrt3−1，面积Q。它包含负y部分，按上一节周期lift计入总高。

该矩形与pair新增底部块的交面积O=(sqrt2+sqrt3−3)(sqrt2−1)。端点次序由sqrt2∈(7/5,3/2)、sqrt3∈(17/10,7/4)、sqrt2+sqrt3>3及相应平方比较严格确定；没有省略接缝或把重叠面积重复算作新收益。

独立展开得

    A123=7/2+2sqrt2−(3/2)sqrt3−sqrt6,
    A123−A^(2)=Q−O=−1/2+4sqrt2−(3/2)sqrt3−sqrt6.

有理根界sqrt2>140/99、sqrt3<26/15、sqrt6<49/20的平方均严格正确。按每项符号代入给

    residual>560/99−111/20=211/1980>1/10.

所以这是完整strip面积的严格缺额，而非只由漏点/小区域推断。它也恰等于完整halo超出全部single/pair union的真实面积；没有pair负Möbius项、重叠或其他面积抵消来掩盖预算。

## 4. U严格见证与全尺度tail

U=(7/16,1/2)×(−1/32,1/32)，面积1/256；固定r=13/8在所有角极值仍捕获三源，平均192/169>1。各单源最小L∞捕获距离下确界为23/16、47/32、3/2，均>1。pair12需半径>23/16>sqrt2；pair23需>3/2>sqrt2；pair13横跨度已迫r≥3/2>sqrt2，故U独立验证纯三源激活。它本身只给集合覆盖失败，面积模板失败另由上节完整积分承担。

r∈[2,13/4)不能跨周期标签，总捕获质量≤12，平均≤3/4。r≥13/4，对3条竖向轨道各用计数≤r/4+1，平均≤3/(4r)+3/r²；该式递减，端点87/169<1。因此包括r=2及所有更大半径均排除，H_R=H。尾界对所有中心成立，不是有限采样。

## 5. 势函数排除范围

单源累计势及每pair的signed Möbius密度累计势，其总charge只分给该pair两个源端点柱。当前柱无任何源，所以全部这些charge为0，与总charge如何signed/asymmetric分摊无关；该边势差必须等于Σpairs Aij−Σsingles Ai。当前m0也为0，任何z均无补偿作用，故上述正residual严格否定此规范模板。

结论不能外推到任意两柱函数、能在空柱收费/重新分配的储备、非累计pair势或三源/整体几何势。源稿明确这些限制，审查通过的亦仅是此限定命题。

## 6. 独立实际精确计算

新脚本 `computation/check_pair_cumulative_triple_review.py`，Python3.15.0标准库Fraction/isqrt，单进程exit_code=0，无随机性、无浮点容差。

算术在Q(sqrt2,sqrt3)的基(1,sqrt2,sqrt3,sqrt6)上精确执行；端点比较用逐步加密的定向isqrt有理包围认证，非零符号确认后才排序。该实现从每个非空源子集bbox和真实平方根门槛生成接收矩形，逐周期lift裁切，以独立x扫描及y区间并求面积。全部7配置/19原始子集的代数结果与上文一致，另独立合并全部purepair矩形求得同一面积4−2sqrt2。

还核对211/1980严格根界、U角极值见证、tail的两个端点、范围增长式的96个有理分支算例。未运行作者脚本或外部大图。

结果 `computation/results/pair_cumulative_triple_review_exact.json`。
复现 `python3 planar-centered-square-weak11/computation/check_pair_cumulative_triple_review.py`。

有限精确计算与全几何证明分开记录；它不证明一般pair函数类不可能，也不给C₂新界。未修改作者稿或第32/42轮冻结文件，未提交、未推送。
