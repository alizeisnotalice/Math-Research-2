# R24 仿射守恒质量族外包络：独立审查

日期：2026-11-15（Asia/Shanghai）。**结论：修订后的盒包络证明、实现及最终分支树通过本次限定范围的独立审查。** 固定族为 204 标签、固定坐标与超周期、参数域 `D=[-1,1]×[0,10]` 的二维仿射守恒质量族。本审查接受下述准确上界作为该族的计算辅助认证结果；不把它提升为一般二维中心正方形弱 `(1,1)` 常数 C₂ 的上界。

终态是 `leaf_limit`：5000 叶、9999 节点，保存目标 `10⁻⁸` 未达到。准确 gap 约 `2.1427153472464046×10⁻⁸`。原始 JSON 的 `scope` 字段仍保留生成时的 `pending independent review`，只是未改写的来源元数据；本审查已对该固定 SHA-256 的终态完成。

本次接续只实际运行单核 `Fraction` 终态树及破坏证书检查，没有重启搜索，没有重算真实 5,081,629 个坐标盒、24,646 个保留盒或真实各节点的矩形并。真实大几何的可信链由完整盒定理接口、源码审查、此前已执行的独立小例及保存运行组成；本审查不是另一套程序对所有真实几何的独立复算。R23 reference 面积也只核对既有证书身份和比值，没有在这里重新计算。

## 1. 固定对象与定理审查

对固定周期 `Px,Py>0`、阈值 `λ>0` 与 cutoff `R`，全部标签坐标固定，质量 `w_i(a,b)=c_i+u_i a+v_i b` 在整个闭矩形 D 严格正，总质量恒为实际值 `M₀>0`。统一尾界在同一 D 上排除所有 `r≥R` 的严格 λ 见证。完整坐标盒定理因而对每个共同参数 θ 给出有限接收矩形并的完整 halo。

每个源盒 B 的质量 `M_B(θ)` 是共同参数的仿射函数。在参数小矩形 C 上，其最大值恰在 C 的四角之一达到。逐盒取此最大值并扩大接收矩形，再取并面积，得到对所有 `θ∈C` 的上界。这一步是合法外包络，不要求不同盒的最大质量由一个真实输入同时达到。分母始终是实际固定总质量 M₀，不能换成各标签或各盒最大质量之和。

半径 `min(R,sqrt(M_B⁺/(4λ)))` 用 dyadic 向外根界替代：整数平方比较核对 `k²B≤A·2^(2bits)<(k+1)²B`，再用 `(k+1)/2^bits` 并 cap 于 R。所有端点为有理数，整数缩放后的并面积可准确求出，未使用浮点根号或浮点面积作证据。

严格开接收矩形使精确跨度等号也为空。实现对**每个 cell**先检查 `λ·span²≥M_B⁺(C)`，再生成 dyadic 矩形；只在全参数域做一次跨度检查不足以保证点 cell 与旧核一致。这是中断前独审发现并已关闭的问题：某点旧核/修后面积为 `17/4`，缺少 cell 检查的初版为 `9/2`。初版多出的面积仍是松上界，但破坏既定逐点一致性和相应收敛解释；以下哈希均指修后版本。

文稿的均值定理常数 `(|β_B|+|γ_B|)/(4sqrt(λm_B))` 正确，cap 为 1-Lipschitz。有限矩形端点连续给出 J 的连续性和紧域达到。固定 bits=40 会留下 dyadic 误差地板；只有同时缩小网格并提高 bits 才能从该证明推出任意小误差。文稿已明确这一条件，没有把本轮固定精度运行写成任意精度完成。

## 2. 真实族合同与统一尾界

中断前已独立从四角逐标签重建全部 204 个仿射系数，不调用新核的匹配或几何函数。已保存的 `real-family-contract.json` 给出：

- 参数域准确为 `a∈[-1,1]、b∈[0,10]`，复制因子为 `2×2`；四角标签质量严格正，最小值 `1/1250`，故整个 D 严格正。
- `Σu_i=Σv_i=0`，实际固定总质量为 `65941670400277227073796193886781311/3125000000000000000000000000000000`。
- 标签短 lift 全部满足固定 `δ=10009397153388108069/409600000000000000000000`，四角最大支配系数 `κ=1001/1000`；仿射性把支配推广到整个 D。
- 按旧基周期和旧基质量代入已审的逐标签输运尾界，准确得到 `κGδ(5)<λ=1`。保存的准确正余量见该 JSON；尾函数对 `r>0` 严格递减，故同时排除所有中心与 `r≥5`。

本次接续把终态 family 元数据逐分数绑定到上述合同，并核对真实输入文件的 SHA-256。该标签与标量核查覆盖真实输入；没有将“全标签核对”表述为“真实全盒重新枚举”。尾界推导的来源接口另见 [R13 独审](supercell-mass-domination-review.md)。

## 3. 已有独立小例与实际执行范围

中断前的 `exact-checks.json` 记录已实际运行的独立核验：直接逐标签、逐周期副本求源盒质量，使用 `Fraction` 端点和另写的 x-slab 区间合并面积，不使用新核的前缀表或原并面积函数。共 98 项 cell 包络比较全部逐分数相等，28 项点 cell 与旧完整 halo 核的 outer 结果相等。覆盖碰撞、wrap、非单位 λ、实际总质量不同于基线质量、非 dyadic cutoff、非标准参数域和退化点/线 cell。

同一记录还含 85 个 pilot cell 的独立重放、单进程四叶分支树、32 项裁剪矩形增量检查、固定 bits 误差地板例、非法仿射/负角点/非守恒/变坐标/尾界失败拒绝，以及零角点模式边界。严格正性是本轮收敛率预算的假设；允许零角点的 kernel 上界模式没有被冒充严格正族模式。

`parallel-driver-check.json` 另记录真实执行的 2 worker、4 leaf 小例，全部七节点回链、相邻二分、叶上界与独立单进程证书一致；越界 reference 在调度前被拒绝。辅助两盒例实际输出准确并面积 9，并明确两盒角点最大值不能同时达到。

本次接续已回读这些脚本、结果和精确来源哈希，没有无故重跑上述小例。并行驱动按 fork 共享同一只读族、完整记录每次相邻中点二分、只在准确 `upper<lower` 时剪枝，最大叶上界提供全域界；wall/CPU 时间只决定停止时刻，不参与数学不等式。

## 4. 最终树与 reference 的新核验

最终原始证书为归档路径 `computation/results/round024/branch_tree.json`，原临时路径为 `/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round024/envelope/bb_pending_audit.json`。本次已只读回查 producer 归档副本，确认其与原临时文件逐字同 hash；归档族输入和两份 R23 reference 文件也均匹配下表。终态固定 SHA-256：

`90b6ad00bff70c4346a91063f7e89588bde8e6e6e1c95e13a9b6699a69781108`。

本次新检查 `affine_mass_final_checkpoint_check.py` 不导入 envelope 或 halo 核，只用标准库 JSON、整数和 Fraction；实际运行 Python 3.15.0、macOS arm64，CPU 约0.265秒、退出码0。核对如下：

- 9999 个节点的 id 连续，除根外每个节点恰出现于一个更早父节点的 children；没有重复、丢失或孤立节点。
- 4999 个 split 均有两个子节点；每对子 cell 正好是父 cell 的某一坐标方向相邻中点二分，回链相符，child upper 不大于 parent upper。完整有限树的归纳证明给出 D 的全覆盖和叶内部不相交；另准确相加叶面积为20。
- 5000 个叶中 4582 active、418 pruned；所有 pruned 满足准确 `upper<lower`，active 满足 `upper≥lower`。保存全局 upper 精确等于全部叶 upper 与 lower 的最大值，唯一最大叶 id 为8911；gap 精确等于 upper−lower。
- `leaf_limit` 与 `leaf_count=max_leaves=5000` 一致；准确 gap 大于 `1/100000000`。没有将失败的目标精度改写成完成，也没有把终态误报成尚在运行。
- reference 参数 `(-1/8,42/5)` 位于 D；其周期、λ和全部204行坐标/质量逐分数等于同一仿射族该点。primary 输入、证书和旧核源码哈希一致，`inner_ratio=λ·inner_area/M₀`，保存 lower 与该准确分数相同。

本次还实际执行8项证书破坏检查：重复子引用、错误父回链、相邻分割错位、错误全局最大值、错误gap、假target状态、越界reference、错误族hash均被拒绝。该检查不承担“篡改某叶几何值”的独立识别；真实逐叶几何值仍来自已审实现。

令 `J* = max_{θ∈D} J(θ)`。准确结论为

```
5787126896633734261086747977203468022712787203324669793218439
---------------------------------------------------------------- ≤ J*
3398607813753623035386245175756010289388060672000000000000000

      5787126969456225479101126541636557987858604811893699578374689
J* ≤ ----------------------------------------------------------------
      3398607813753623035386245175756010289388060672000000000000000
```

上界约 `1.7027934044159632`，reference 下界约 `1.7027933829888098`；小数仅展示，认证使用上述分数。准确 gap 为

`149140462014493447299958968248618634462349373/6960348802567419976471030119948309072666748256256000`。

## 5. 来源身份、交付与复现

| 来源/产物 | SHA-256 |
| --- | --- |
| `proofs/affine-mass-rectangle-enclosure.md` | `b67ae00033b72fc991a8a679510c8de2ea0f3ccef353161d712e93d7820152fd` |
| `computation/check_affine_mass_enclosure.py` | `554f8694a7cc188ef9c62807d1b2c1920777c353d9f05bdbbae2a6be45f36b0a` |
| `computation/mass_family_envelope.py` | `19db8ef78e753b3d476b66aef4ce435c54d9f8c09c0e22c4664eb4f954bfb5d7` |
| `computation/parallel_branch_mass_family.py` | `a9e7cc66bf2ac79494a425064aa5432bd4ea52b8ef57f141da53016892623096` |
| `computation/periodic_halo.py` | `b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33` |
| `computation/inputs/r024_common_family.json` | `cc43a36cd1fa0eeb486fe251fb6b47eca539abcae9ca402c596f14feea9f2322` |
| `computation/inputs/r023_final_1051.json` | `4fb1461289fd2b1ef52ac7d30c3d5cd0bd5b72c102f7ccdbeb01be6706baf0ee` |
| `computation/results/round023/primary_final.json` | `909c10068b2122b4c514ce1ff91f2886c7c5aae58c02e3cd050fe73339600a99` |
| `computation/affine_mass_final_checkpoint_check.py` | `c4c15acd59f34fc55cade9727f36a9aa01c903c3dbc649a8636660155a59033e` |
| `computation/results/affine-mass-envelope-audit/final-checkpoint-check.json` | `f1ac339db17db3e2c8aadc4a6c020a1bed67f3141ec44dd4e703eb8d8d6c13b5` |
| `computation/affine_mass_checkpoint_rejection_check.py` | `71a19c68c398aad44b3a1ed3f4ea6d8503bd344d08dddba1f25b6c4b4b2dd3e9` |

前阶段独审输入、脚本hash、具体准确结果与复现命令保存在 `computation/results/affine-mass-envelope-audit/{exact-checks,parallel-driver-check,real-family-contract}.json`；新终态与拒绝检查输出分别为 `final-checkpoint-check.json`、`final-rejection-check.json`。既有 `live-checkpoint-check.json` 只检查当时3397叶快照，保留其历史身份，不再拿它代替终态。无随机算法、容差或新数学库依赖。

终态检查器的每个 `--参数` 路径和实际解释器均已保存在新输出的 `actual_command`。归档后可把其 `--snapshot/--family/--reference/--certificate` 换成上表稳定路径，其余参数保持指向同hash的证明、实现和前阶段独审JSON，再运行检查器。拒绝检查可复现为：

```sh
python3 -B computation/affine_mass_checkpoint_rejection_check.py \
  computation/results/affine-mass-envelope-audit/final-checkpoint-check.json \
  /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/affine-mass-envelope-audit/final-sensitivity \
  computation/results/affine-mass-envelope-audit/final-rejection-check.json
```

上例从课题目录执行；拒绝检查按前一次输出的命令定位来源，若实现已从独立工作区移到归档目录，应先按新的实际路径重新生成终态检查输出。

所有新增审查文件只写入独立 `pcw-kernel-audit` 工作区；producer、main 和旧 T 未被本审查修改。本附审不单独提交/推送，由父任务归档 R24。尚未解决的是该族最大值更精确定位与目标 gap，以及远超本族范围的一般 C₂ 精确值；本轮没有核验备份完成或一般输入上界。
