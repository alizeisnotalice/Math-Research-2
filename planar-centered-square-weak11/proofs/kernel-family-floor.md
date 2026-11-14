# 第9轮：固定1182变量核族的严格成本下限

日期：2026-11-15。结论为核族成本约束，不是 C₂ 的下界或新上界。

## 1. 准确命题

冻结 α=47/40、R=33/20、h=1/32。令 I 为
`computation/results/kernel-interface/floor-neumann.json` 的 `variables.classes`
列出的1181个对称样条类，逐项对应旧候选
`88af41536e1117d953bbc4594ee99f1319b79b75d98d740c5d3c82caaf1f306f`
的全部系数键。令 O(i,j)={(±i,±j),(±j,±i)}去重，定义

    Φ₀(u,v)=((|u|+|v|)^−α−R^−α)_+,
    Φ_ij(u,v)=Σ_(a,b)∈O(i,j) B(u/h−a)B(v/h−b),
    K_θ=cΦ₀+Σ_(i,j)∈I a_ijΦ_ij,
    c≥0,  a_ij∈R.

样条系数允许负；c和全部1181类系数均可变化，没有删除变量或固定旧系数值。不得把结论扩展到增加新样条类、改变α/R/h、改变奇异型或使用另一生成元的更大核族。

**命题。** 此族内任意满足核比较完整连续域条件 C1、C2、C4 的核，其半积分 J 满足

    J=(1/2)∫K_θ ≥ L > 3.589342358790,                    (1)

其中主认证可取准确有理数

    L=1408166696828769191441564314799163448294543331165700878499 /
      392318858461667547739736838950479151006397215279002157056.

该L略大于3.58934235879042629824443857909495788423660003729。
独立程序所得下端略宽，仍严格大于(1)中的短十进制有理数。

与本轮重新计算的当前旧核 J₄ 比较，有准确分数检查

    J₄−L < 0.001523.                                     (2)

因此本固定族内继续任意调整全部系数，从当前核出发的成本改进总量不足0.001523。这里旧核的**连续可行性**仍来自历史完整C4认证，本轮重算其积分；(1)的族成本下限本轮新认证且已独立复核。

这并不证明真实C₂低于该成本地板。约1.7的已知构造是C₂下界，不能与(1)相减后宣称已证明核方法与真实常数存在差距。

## 2. 完整变量的有限必要约束

在准确点上取下列必要约束：K≥0；在菱形 D={|u|+|v|≤1}内K≥1；A_αK≥0。C1/C2点避开原点，C4点另避开坐标轴与u+v=R。

有限样条核在原点外连续，故几乎处处C1/C2也强制这些点的不等式，包括菱形边界上的极限。有限样条生成元在所选非轴、非截断界面处连续；分布C4及其点态表示因此强制相应点值非负。点约束是真正必要条件，而非只验证Bernstein系数的更强充分条件。

写1182维参数θ=(c,a_ij)，行矩阵Gθ≥b。积分成本向量为

    m₀=αR^(2−α)/(2−α),
    m_ij=(h²/2)|O(i,j)|,
    J=m^Tθ.

全部变量列都在矩阵中。类轨道数为1个单点、99个四点轨道、1081个八点轨道；不同轨道互不相交。

探索先用3368行必要约束，浮点最优值仅约1.1524，没有改变路线判断。随后使用10657行：5328个K约束及5329个A约束，含完整knot点、中点、轴旁点及旧失败点。浮点LP约53.0秒给出候选目标3.589342357805082，选出1182个正dual行；浮点值本身没有认证作用。

最终仅保留这1182个准确活动约束作为证书输入：15个K≥1、465个K≥0、702个A_αK≥0。逐项点、类型、右端b以及全部变量类均保存在 `floor-neumann.json`。全10657行菜单不是证明依赖，可按脚本配方重建；庞大的浮点矩阵仅存CodexTemp。

## 3. 正确的对偶方向：自由系数列必须是等式

令A=G_active^T为1182×1182矩阵。如果存在y≥0且

    Ay=m,

则任何满足必要约束的θ都有

    m^Tθ=y^TG_activeθ≥y^Tb.                              (3)

这里没有对θ的符号作额外假定；尤其适用于任意正负样条系数。不能把列等式随意换成G^Ty≤m：对于负的a_ij，乘法方向会反转。所有成本列正也不能修复这个错误。把近似有理dual统一缩小，会破坏自由变量列的等式，故本稿没有采用这种伪证书。

探索的有理向量y₀只是精确等式解的近似。下面用只含矩阵乘法的残差定理证明**存在真正的非负等式解**，从而严格调用(3)。

## 4. Neumann证书，不依赖求逆程序正确性

保存准确有理向量y₀及一个二进制浮点矩阵R；R的每个binary64条目均按其准确dyadic有理值解释，不把R当作可信逆矩阵。按准确点闭式，用192bit Arb向外包围实际矩阵A和成本m，实际认证

    q≥||RA−I||∞,  q<1,
    η≥||R(m−Ay₀)||∞,
    δ≥η/(1−q),
    min_j(y₀)_j>δ.                                      (4)

RA可逆，因此方阵R、A均可逆。唯一真解y*=A^−1m满足

    ||y*−y₀||∞≤δ,
    y*>0,
    b^Ty*≥b^Ty₀−||b||₁δ=b^Ty₀−15δ.                    (5)

将(5)代入(3)即证明(1)。所有自由变量列连同c列都参与RA乘法及残差，未做遗漏列定价。

主程序实际数值包围：

    q ≤ 0.0007500577878818873702722876426892260468,
    η < 5.245×10^−47,
    δ < 5.249×10^−47,
    min y₀ > 1.3204×10^−9.

这些十进制仅便读；JSON保存准确向外端点，最终比较用有理数。该认证8.65秒；从保存的最小证书重新执行、不读取浮点LP菜单或旧候选系数，另一次实际复现约5.55秒。

先前的 `arb_mat.solve(precond)` 严谨解包围保存在 `floor-certified.json`，总计17.60秒；它用来发现很准确的y₀，但不是(4)–(5)验证的逻辑依赖。即便该求解器输出任意y₀和R，只有(4)成立才能通过最终验证。

## 5. 独立实现复核

另一个工作区的独立程序没有导入本程序或vendor `hl_core.py`，未调用solve/inverse，而是直接按轨道展开、独立重写内外奇异闭式并计算完整1182矩阵。它核对全部类、活动点必要性和(4)，实际约6.66秒得到

    q=0.0007500577878818873702722876426892260467…,
    η<5.271×10^−47,
    δ<5.275×10^−47,
    b^Ty₀−15δ>3.58934235879042629824443857909495788423660003728.

独立输出位于独立工作区 `pcw-kernel-audit` 的
`computation/results/kernel-interface-review/floor-neumann-independent.json`，
脚本为 `computation/kernel_floor_neumann_review.py`。父代理负责白名单整合；本稿没有改写独立代理的材料。其证明方向与本稿相同，但矩阵及闭式为独立实现。

## 6. 必要预条件输入的归档与重现

原R为1182×1182、little-endian binary64、row-major，共11176992字节，SHA256为

    2c6e75a224177c81875f78dce3d5d314f25e1da03a1e6eb518c99821f87f2dfc.

普通gzip大小10810771字节、xz大小10596156字节，均略超10MiB。采用可逆8字节plane shuffle后gzip(level9,mtime=0)，保存为

    computation/results/kernel-interface/floor-preconditioner.byteshuffle8.gz

大小9822059字节，小于10MiB，SHA256为

    6765a557999b826e7910f14cfd9be485b1ce3864544f2a10869a66d61952e937.

实际逐字节往返核对成功。`floor-preconditioner-storage.json`记录编码、大小、两项哈希。这里保存R是为了独立复验取得同一必要输入，不把它作为可信逆矩阵，也不归档探索浮点矩阵缓存。

以下均从课题目录运行；`PY_ARB`指已有含python-flint0.9.0的Python环境。当前实际环境是Python3.12.13。本地依赖读源码的vendor已随代码按Apache-2.0保留；最小验证不再依赖旧T。

```sh
python3 computation/unpack_kernel_preconditioner.py \
  computation/results/kernel-interface/floor-preconditioner.byteshuffle8.gz \
  /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-interface/floor-preconditioner.bin

"$PY_ARB" -B computation/kernel_family_floor.py certify --neumann \
  --root computation/kernel_vendor \
  --input computation/results/kernel-interface/floor-neumann.json \
  --preconditioner /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-interface/floor-preconditioner.bin \
  --out /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-interface/replayed-floor.json
```

`unpack_kernel_preconditioner.py`只需标准库，先核压缩文件哈希，再恢复8个byte planes并核原始R哈希；独立checker无需改变已审代码。

保存的准确点如何生成探索的A_float和R，亦有可重现命令。`PY_SCIPY`指已有NumPy2.5.4/SciPy1.18.1环境：

```sh
"$PY_ARB" -B computation/kernel_family_floor.py prepare --step 1 \
  --root computation/kernel_vendor \
  --certificate computation/results/kernel-interface/floor-neumann.json \
  --out /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-interface/rebuilt-menu.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$PY_SCIPY" -B \
  computation/kernel_family_floor.py preconditioner \
  --input /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-interface/rebuilt-menu.json \
  --lp computation/results/kernel-interface/floor-lp-dense.json \
  --out /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-interface/regenerated-R.bin
```

`prepare`给准确点及Arb singular列；`make_matrix`给探索浮点矩阵，按保存的活动行构造A_float，再以SciPy求R。平台变化可能改变R的比特，故新的R仍必须重新通过(4)。归档的原R则由上述固定哈希锁定。

以上两条生成命令已从随包最小输入实际重跑：菜单重建约2.94秒，R生成约1.31秒；新R的11176992字节及SHA256与归档原R完全一致。解码归档R后的最小验证也实际重跑约5.61秒通过，摘要见 `floor-bundle-replay-summary.json`。这些重现不依赖旧T文件或LP求解器的最优性判断。

## 7. 研究决策与未解范围

本轮把“继续调固定1181类系数”从开放盲搜压缩为剩余成本空间不足0.001523的受限问题。下一次希望有显著一般上界改进，应明确改变冻结条件，例如同步改变奇异幂与生成元阶数、改变R/h或扩展样条支撑/基函数；每次仍须新的全连续域C1–C4认证。

本轮没有找到新可行核、没有重放历史完整C4、没有算出该核族准确最小值，也没有确定真正C₂。低成本浮点LP候选不是一般上界；固定核族的严格成本地板也不是C₂下界。
