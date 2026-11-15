# R49交叉独审：排序障碍与固定八块指定序上证书

2026-11-15，Asia/Shanghai。审读者工作区pcw-sol-cover-structure。结论：所列限定命题通过独立数学审读；唯一指定八块顺序的严格outer面积已独立准确重算。没有一般K3/K3.5新结论，未改R34/R40、作者来源或生产输出，未提交/推送。

## 绑定的作者最终来源

来源根为 `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-cover-counterexample/planar-centered-square-weak11`。以下均实际读取并按字节核对哈希，未仅信manifest中的PASS。

| 文件 | bytes | SHA256 |
|---|---:|---|
| `proofs/post-r40-sorting-obstructions.md` | 11566 | `47f0ac8ff9675a02458ef313654f7f44b0153fb7fd55e4cd1d340313f1e08beb` |
| `proofs/quarter-cell-prefix-group-single-test.md` | 3984 | `8faa60748b4215b8bd8662c3047e31da7ddaf106bd9801d358e53fec17f1f510` |
| `computation/post_r40_sorting_obstruction_check.py` | 5215 | `e533adf57f104bcf822df6fff2eafab81874d92711b7fcd8efcb433123b20291` |
| `computation/results/post-r40-sorting/exact-obstructions.json` | 4557 | `cf18aac220582195958c5dce0b20a153ada2a5e228ec8009430d3f44667a3ff8` |
| `computation/results/post-r40-sorting/HANDOFF-MANIFEST.json` | 1625 | `4047e8c0770f5fa53bf28e32451d3b3b3ae5cd70a5a1b8169e3b99de1f1dfac8` |
| `computation/quarter_cell_specified_order_upper.py` | 2904 | `bfcb7fe72e550cf3f9eb7d04e3fc9994ef30c2b58ad91390b3279a03108587ee` |
| `computation/results/quarter-cell-prefix-group-proxy/quarter-cell-single-input.json` | 2611 | `381cc5adc4516b2fd86e6ec4fca8d7b02114f8880c9f14311257402533350299` |
| `computation/results/quarter-cell-prefix-group-proxy/specified-order-outer-certificate.json` | 4290 | `67bed2dbb6eaadb528d3d308c5843c914c149f0721652c97e4bfbc6413d3d23a` |
| `computation/results/quarter-cell-prefix-group-proxy/HANDOFF-MANIFEST.json` | 2331 | `763e78584d54cd20ce6a809831d95d481536407a9751a4cc9750dd9c978e5dbc` |

## post-R40排序障碍：数学核验范围

1. 对任意非零线性functional，八向量分为四条不同无向直线的反向对，kernel最多含其中一对；至少3个严格正projection轻点必晚于重源，与tie如何处理无关。每个后轻点+重源eligible，因此四个unit squares被proxy包含。不同grid中心∞距至少99/100，任意pair交面积≤1/100，Bonferroni给197/50；除总质量1251/1250得4925/1251。这个论证连方向按输入选择也覆盖，但只针对单投影排序。

2. 重源最后的精确自由最优10011/10000核准：全部轻点组质量小于最小pair跨度平方，所以轻阶段全singleton；任何4轻点bbox maxspan=2a，因而重源可用恰3轻点。半径和严格小于间距，得互不交；重源非最后则至少两unit squares union≥199/100，严格较坏，故这里确为全序最优，而非只找到一个好序。独立脚本实际检查70个四轻点bbox和准确比值。

3. 100中心原子个体质量升序障碍核准，包含任意ties。中心组diam≤99/100000，末中心proxy含unit square；每个后perimeter点共同捕获全中心质量的跨度≤99099/100000，其平方<51/50。八外围unit squares只留中央开hole(−49/100,49/100)²，任何允许中心unit square都覆盖该hole，故union至少(149/50)²，除29/25得到22201/2900。这里unit squares是较大proxy内的保守子集，不把proxy错当面积恰1。

4. perimeter先、center后的好序上界核准：前八proxy total4/25；任何中心有效组至多含3外围源，所以所有中心proxy area≤53/50，其中心x-span<1/1000，union装在宽√(53/50)+1/1000、高√(53/50)的矩形内。用√(53/50)<103/100得到总<122103/100000、ratio<122103/116000。不是精确自由最优声明。

5. 三点角簇交错例严格阈值核准：p0,m1；q1,m1；r4,m14。q eligible，r和q+r均不eligible（后一是16=16，不能用≤）。但r→q使q proxy area15，包含p proxy area2，外损13。它只否定任意interleaving下向C单独收费，不否定存在性联合预算。

6. 组件独立及条件separator预算、rank≤k按新源质量付款、proper rank2根大边预算均核准；它们没有假设P=M。六点环pair graph C6、删点后path5连通但整体非eligible、proper rank2由全bbox比较成立；它因此排除普遍小separator，却由rank-root预算8/3处理。独立脚本检查全部15 pairs、41个proper≥3子集、6个删点bbox；连通性由明确cycle删点为path的解析事实承担。

## 固定八块容器：边界、任意块内源与作用域

作者文中的闭区间[0,1/4]、[1/4,3/4]、[3/4,1]应理解为cell的闭bbox容器，而非字面互不重叠的measure partition。对块内任意源的类，必须预先给源质量唯一归属，或写为8份各质量1/8的labelled子测度ρ_i，每份支撑在相应闭bbox；共享cut的质量不能重复计算。闭bbox只会扩大接收容器，独立几何证书不受边界源位置影响。

在这个严格解释下，作者“覆盖块内任意源”成立，不局限于某个原子网格。对有限原子，按块顺序插入，rank k块内任意prefix的总质量≤k/8，所以每个point proxy半径≤√(k/8)/2且中心在该块，故U⊂G。甚至对任意这样的有限正块测度也可直接验证真实halo包含：给实际strict witness取捕获正质量的最后block rank k，则捕获质量≤k/8，故r<√(k/8)/2；该方形与B_k有质量交集，中心在B_k+[−r,r]²⊂R_k。于是Hρ⊂G。这里没有声称G等于U或真实halo。

这个“任意块内源”只对应固定八块总质量向量(1/8,…,1/8)、中央质量0及固定几何。不能用它推出所有c>z、所有九格质量或所有宽大bbox输入的3.5预算。

## 唯一指定序严格证书的独立复算

绑定序(5,6,7,3,0,4,2,1)，rank by cell=(5,8,7,4,6,1,2,3)。逐个独立验证8个upper radius numerator H满足(H−1)²<k·2^59≤H²，因而H/2^32确为向上ceil根界。全部cell端点乘D=4·2^32后为整数；重新构造8个outer rectangles，与来源证书逐整数相同。

作者使用2D endpoint-cell覆盖计面积；本审查没有调用/导入该helper，而以独立vertical sweep在每个x slab合并y intervals，重算准确整数并面积

943119706436290073136 / 295147905179352825856
=58944981652268129571 / 18446744073709551616 < 7/2。

严格margin为5618622605715301085/18446744073709551616，与来源一致。证明方向正确：outer rectangles包含真实根式容器，故其并面积是上界，足以证明存在这个序低于3.5；没有认证10080浮点探针的精确最优，也没有要求证它最优。

实际运行：Python3.15.0标准库Fraction/整数，单进程，无随机数或浮点决策，exit_code=0。耗时0.004451208020327613秒。脚本为`computation/review_r049_counter_and_quarter.py`，实际输出`computation/results/r049-counter-quarter-independent.json`。输入只读取作者冻结文件；未向作者目录写入。

未重跑：10080个浮点排序探针、生产2D endpoint-cell union、16-sector原作者checker、任何大搜或MATLAB。本文全方向/任意tie声明来自上述解析审读，不把有限数字检查替代全量词证明。

## 与本批四边账本的边界一致性

本作者自己的结构稿末尾追加固定ambient[0,L]²的四开edge条带与闭center[a,L−a]²，逐点h+v=1+hv−(1−h)(1−v)，得到S=M+c−z。此为边界口径核对，不冒称对本人稿的独立审查；其正式独审由counterexample协作者承担。该账本不能用本节闭bbox容器的overlap重复计质量，也不能混同实际bbox宽度不足2a时会互相重叠的opposite strips。

## 最终状态

两份作者证明的限定结论与指定序outer证书通过。本review不覆盖一般九格mass simplex、一般free-order K3/K3.5或精确C₂。后续九格与跨大bbox组件工作须使用新R53材料；本R49交叉审查和已锁定来源停止编辑，供root归档。

