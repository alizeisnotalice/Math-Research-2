# 第49轮结构稿交叉独审：通过

日期：2026-11-15（Asia/Shanghai）。审查方 `sol_cover_counterexample`，作者 `sol_cover_structure`。**结论：最终锁定版本的限定命题与证明成立，无数学修改要求。** 本审查没有证明一般 minπ|Uπ|≤3M、一般3.5M或精确C₂。

## 最终版本绑定

作者工作区为 `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-cover-structure/planar-centered-square-weak11`。本次完整读取、独立推导并在检查前后复核以下最终三件 SHA-256：

| 文件 | SHA-256 |
| --- | --- |
| `proofs/sol-effective-separator-recursion.md` | `a913f6eb302be63208214a92ae858be7e13a95d4561b1eaf5b228a6149bb95b5` |
| `computation/sol_effective_separator_check.py` | `260908fc99bb3415aa0757f634907bffe3acaa5818e2eda58e0217dea1dc08c0` |
| `computation/results/sol-effective-separator-exact.json` | `6350998ebd1c149ae3563750b5f9071f53332bf2531738cc6b705177b00168e8` |

没有把此前中间版本或旧“通过”标签当成本次审核。本次未执行作者脚本。作者结果JSON中的 `proof_sha256=fb7cf579…` 对应追加最后纯解析附注之前的实际运行，不是最终正文的重跑记录。最终附注由本次直接审读与下述很小的边界算术检查承担；保留历史运行哈希不造成证明范围混淆。

## 逐项数学审核

**严格eligible与角簇。** d(B)²<m(B)恰对应闭共同半径见证可具有严格质量余量。固定p为同一闭象限角时，d(p,S∪T)=max(d(p,S),d(p,T))，结合非负质量，确实保持严格eligible；有限并给唯一最大簇C。对非空外D，C∪D不eligible且d(p,C)²<m(C)+b，严格推出d(p,D)²≥m(C)+b+m(D)。等号必须保留为≥。三点p0/b1、q1/m1、r4/m14例正确，C≠旧eligible组件，任意interleaving的内部外损13不能向C质量1独立收费。

**eligible组件与完整源盒。** 每个prefix winning标签组也是完整旧输入的eligible超边，因而不跨该超图分量；反向候选亦保持，故分量内proxy逐步完全相同。真实halo见证也只捕获一个分量；Hρ=并H(C_l)成立。补入bbox内全部源不扩大bbox却增加质量，因此完整源盒可生成同一组件；这不要求相交eligible边取并仍eligible。

**删点及mask递推。** 最后Q_p固定，旧源可按eligible组件独立排序。式(3)的≤正确，因其只丢弃组件在Q_p外的接收交叠。n点质量升序的proxy面积和≤(n+1)m/2，给出≤3源与≤4源组件的两个受限结论。一般B(A;W)的式(5)是F(A;W)的**上界**：对每个p构造孩子顺序，取并后以面积次可加性支付；对孩子递归归纳即可。不同孩子接收集可能重叠，所以不能把该递推写成F的等价定义或已经求值的一般常数。遮罩须累加W∪Q_p，源标签在父删除点和孩子组件之间不重复。

**低rank、proper-rank根、weighted低cap。** mass非降序使winning边每个原子质量≤当前新原子质量，故rank≤k确有T_j≤km_j。最大质量点最后时，其前winning边全部proper，给式(9)。当前子问题若可删点T_p≤3m_p，配合不交组件和成本≤自身质量的leaf，确能以删除点或leaf内标签分别付款；稿件未把这个充分条件冒充总可选定理。六点环既无≤4组件separator，又满足proper-rank2，故8/3结算正确；“大组件存在”本身不是3M失败证据。

**dense leaves。** 若D≤(√M−√R)/2，旧proxy半径≤√R/2且中心距p≤D，故全部落入最终S_M(p)；全源bbox严格eligible保证最后T_p=M。即使条件取等号，开proxy仍保持所需严格内部距离。更强的4D_bbox√M<m_min规则同时保证所有完整prefix eligible，以及每步半径增长≥m_min/(4√M)>D_bbox，从而任意顺序nested。带遮罩时前者是构造上界，后者因所有顺序都退化为某个最后源方形，才有相应min_p等式；二者区别正确。

**单条带停止规则。** 左开条带先插入时半径≤√R/2；其余点在cut上或右侧，半径≤L/2。条件a+√R/2≤L/2正好把两阶段左边都放在x_min+a−L/2以内，其余三边由bbox和总质量直接控制。条件不要求W,H≤L；只有导出3.5M的特例才加这个假设。cut原子归后阶段，因此等号R=M/4与a=L/4允许。完整源bbox严格eligible才得到严格<3.5M。三点样例容器面积17127/5000与实际proxy3/2分别正确，没有把容器面积称作实际halo面积。

**真实bbox birth接口。** 含p的全部源子集接收矩形并R_p，加上删p旧组件的真实halo，确给H(A)的精确并式；每个非空矩形来自eligible标签且包含p的距离约束，所以R_p⊂Q_p。该公式没有覆盖原proxy U的承诺：用R_p替换Q_p只能服务实际H的账本，不能据此缩小原U后仍宣称覆盖它。稿件已经清楚保留此范围。

**四边与四角账本。** 最后附注的ambient是固定[0,L]²，而不是可能窄到使相对条带重叠的实际bbox。开边条带与闭center定义完全匹配；h,v逐点为0或1，含cut原子时仍有h+v=1+hv−(1−h)(1−v)，故S=M+c−z精确。角质量在四条带中被计两次，不能当四份互不重复储备。K∈[7/2,4)时a=(4−K)L/2、R≤(K−3)²M确给a+√R/2≤L/2；直接按ambient端点重做条带容器，宽≤2L−a、高≤2L，面积≤KM。这个应用只用源被ambient包含，不要求ambient恰等于实际bbox。全部四条带失败才推出c−z>[4(K−3)²−1]M；K=7/2仅为c>z，不是一般闭合证明。

## 本次实际计算与复用界限

我方实际运行 `computation/review_r049_structure_final.py`，Python3.15.0，Fraction，单进程，退出码0。最终运行耗时约0.018709秒。它仅核对：最终三件哈希及作者列出的四个既有main依赖哈希；a=0、1/4、1/2下43个含cut/角/中心边界点的逐点恒等式；4个有理K的标量代数；一个条带容器标量。有限点检查不替代上面的全范围解析审读。

已复用、未在此重新执行的是我方 `computation/results/post-r40-sorting/exact-obstructions.json`，SHA-256为 `cf18aac220582195958c5dce0b20a153ada2a5e228ec8009430d3f44667a3ff8`，它记录之前实际运行的六点环/41个proper子集、投影与interleaving小例。作者脚本只是完整阅读，未重跑其subset或proxy检查；源JSON作为作者原实际运行记录读取。

作者脚本source索引的四个依赖在当前main均已存在，且本次实际哈希与作者输出完全一致：`sol-effective-envelope.md`、`sol-effective-ordering-draft.md`、`general-dual-interface.md`、`mass-capacity-equality-attempt.md`。本审核没有因此自动重证它们所有历史结论。

复现本次最小审查检查（整合到main后，默认author-root即当前topic）：

```sh
python3 -B computation/review_r049_structure_final.py
```

在独立clone使用 `--author-root` 指定锁定作者topic路径；实际本轮另用 `--main-topic-root` 核对main既有依赖。实际命令及路径保存在输出JSON。

## 最后白名单与未解问题

作者三件、我方本review/检查脚本/实际结果/小manifest为新增白名单；具体SHA与大小见 `computation/results/r049-structure-review-manifest.json`。四个已在main的来源文件不需要重复复制；复用的post-r40结果应随该作者包保留。本次没有修改作者源、main/T或R34/R40及已冻结R49正文，未提交、未推送。

仍开放的是一般高rank重叠结构的可传递源质量/接收遮罩势、所有输入的一般3M或3.5M，以及真正C₂匹配界。正式审核通过的对象是本稿全部明示假设下的命题和有限接口，不是这些尚未证明的总目标。
