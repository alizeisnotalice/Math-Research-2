# R39 完整单轴slab枚举器独立审查

结论：PASS（完整枚举器及所列真实输出核验范围）。绑定 `computation/position_slab_enumerator.py` SHA256 `acede2c710655d0c179e7bea3ee9483271caaf37968227e7e9fce6442907056b`；真实 `slab.json` SHA256 `c0d262c39a01db8d50cd2fee1a00c7867aa8e08417c4fc85522103b70ce13b25`。未改作者代码，未重审冻结理论，未回放44,045,029盒的大循环。端点排序、扫线β和最终δ组合由envelope另审，本报告不单独认证β>0。

## 源码与分类

已独立读完整实现。周期副本I按初始[s₀−H,s₀+H]对窗口的可能相交范围一次建成，之后只缩H、不丢冗余副本；零质量标签删去合法。源x按(x,velocity)合并持久类，不同速度在中心碰撞时抛WallAtCenter。全部source_order、span=2R、源可达F±R的非持久事件在任何span/reach预筛选之前处理，H取半事件距离的最小值。故按中心筛掉的坐标对在整段保持失败，bbox成员和质量固定。全部有序端点对含退化坐标对，整数prefix使用闭盒下标，碰撞质量相加正确。

质量Q的分母清除已逐式核对。记坐标公分母S、质量公分母T、λ=ln/ld、D为整数x跨度、M为整数盒质量，则Q=ln·D²·T−M·ld·S²，物理q=Q/(ld·S²·T)。代码给出的排除半径准确等于q/(2L)，L=λ|u|(2|D|/S+|u|H)。clip分支采用有效半径平方r²=K/(B·S²)、B=4Tln，Q_clip=dist²·B−K，物理平方余量Q_clip/(B·S²)，其半径也准确等于q_clip/(2L_clip)。未混用λ、质量或坐标分母。所有计算均用Python整数/Fraction，无浮点安全余量。

分类穷尽：零质量；y跨度/固定y裁剪恒空；x质量跨度恒空或有正余量；cap跨度恒空；左右x裁剪恒空或有平方余量；余下才生成真实活盒。根号等号且速度非零会fail closed，持久等号恒空；q>0分支以有理余量保证整段不激活。左右clip平方判据仅在相应距离非负时使用，避免平方产生反向判断。所有访问bbox只落一个拒绝类或活盒，计数恒等式检查该穷尽性。活盒只表示中心当前非空；其后是否出生/死亡、裁剪速度是否切换仍由独立端点全预序证书限制，源码已明确此边界。

## 本次实际小例

Python3.15.0，单进程；作者自检已实际重新执行，exit_code=0，输出与原selftest.json逐字节相同，SHA `bb6b56a09e11a0a369dfbc9f8dd4ee85f4bb8563a7517bd40c253aaa238513bd`。包括9/16两原子、source/死亡/span墙、clip及周期目录边界，以及seed39015的1260次直接接收距离—累积质量oracle。

另以自身脚本作独立检验，实际0.164842秒、exit_code=0、seed390151：8个有限及3个周期有理输入，含λ=2/3、5/7、复杂坐标/质量公分母、零质量标签及未wrap周期位置。独立重建11份完整副本目录，逐点闭盒相加而不调用生产prefix；3834次非零bbox半径/裁剪判据检查，中心真实活盒集合完全一致，在±δ_empty处没有被拒盒新激活。每类正余量见证核对radius=min(H,q/(2L))。另7例准确验证source_order、x_mass_span、span_2R、reach_left/right、x_clip_left/right均fail closed。有限检查辅助完整源码审读，不作为任意输入定理的数值证明。

## R39 真实输出

本次独立只读复核耗时0.993998秒，exit_code=0。输入为R29 candidate_0298，文件SHA `ca9f8b0ea3d78036151ed3e5c37ffd558c395fc2de66c273f561f4ef86630f59`，canonical SHA `7430136b07f641f39e4e26cd12fee6a6fbb5a5ca317a688f4848d90ea03666fa`。速度文件SHA `09b06a3ef8f08a9249a712462d7f5f7f9035b30fc7598c7471acf51c7ae1c2ff`，204标签中仅130/131速度1，其余0。准确s₀=`5588494330693/500000000000000000`，H=δ_empty=`328734960629/500000000000000000`，R=10；目录和参数均匹配。

以原始输入与初始H重新生成全部3953副本、225个x类、62个y类及哈希；独立检查50850条source/span/reach事件，全部支持实际H。重建全部坐标对预筛选：23821个x对、1849个y对，乘积44,045,029，预筛选前49,655,025盒；和输出访问计数及活盒+各拒绝类总数恒等式一致。

用与作者相反方向的y-major整数prefix重建质量，逐项核验全部25,466活盒的bbox成员质量、半径平方、横速度和四个原始端点，并独立检查中心真实裁剪矩形为正；无重复bbox。全部9类拒绝汇总见证的盒质量、端点、速度及实际拒绝条件准确；3类可动余量的q、L和radius均准确，均支持δ_empty=H，最小见证为source_slab。普通周期尾界从当前实际质量和周期重算为

`23820787253964867284592995960140706111190198205884874124835164741/24052915052186143128797603913750000000000000000000000000000000000<1`，沿整段位置变化统一有效。

这些真实检查未逐盒独立重算全部44M拒绝分类；其全覆盖/分类正确性由本次源码审读、全候选对计数、实际记录和上述独立有限oracle共同支持。不能把摘要哈希或count恒等式单独说成每个拒绝盒已独立重验。原真实枚举记录24.730344秒是作者本轮运行，非本独审重算。

## 当前最小端点半径的解释边界

另以准确Fraction实际核对`slope.json`（SHA `479524aff5a851e23f256198c818a77678ca4de7698d957e9eec8f104c7ed597`）的limiting_pair：两个横端点分别唯一来自零基box25244的左端点（bbox[137,138,25,25]）和box25243的右端点（bbox[137,137,36,36]）。两盒固定y中心分别为 `35136509207263900091257231773/64000000000000000000000000000<3/5` 和 `316077158213425564705302455727/64000000000000000000000000000>24/5`，两半径平方均严格小于1/100。故前者y支持严格在7/10以下，后者严格在47/10以上，两支持间距>4。

该最小δ的限制来自两个不相交水平带的raw横端点间隙；它是保持全局raw端点序的保守限制，不能据此声称已定位真实面积饱和墙或最近有贡献事件。此补充不认证完整扫线β，也不排除同一未来参数处发生其他事件。

补充脚本 `computation/review_position_slab_limiter.py`，SHA `4239278d5d5ed2a2a3a9f813b55c8597d0303ce9604b3925ff029799f8f3d73c`；结果 `/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/r039-slope/independent-limiter-note.json`，SHA `2f81819c4019a8acf989a6833cea2c93a3284c27f672ac6004224ab6bde2226a`。该次只读匹配与平方比较，未作大枚举。

## 文件与复现

以下均仅标准库依赖，保存实际版本、命令、种子/参数及输出范围。

- 独立小例脚本：`pcw-r001-audit/planar-centered-square-weak11/computation/review_position_slab_enumerator.py`，SHA `cfa762e57065eccb60e1bb163e604eaf571368906069c178df4520cc596ea2f3`；结果 `/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/position-slab-impl/independent-enumerator-review.json`，SHA `ba7b1dfcf5c300146c8a9b0f4bac5ce8689c013db1875c87f47ea85ada0cbe5c`。
- 真实只读检查脚本：同workspace `computation/review_position_slab_real.py`，SHA `650ff0eff8d5265cfc9c0e7d37e9d06847969b8976abd3d4567a378dc279a599`；结果 `/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/r039-slope/independent-enumerator-real-review.json`，SHA `0a639c5977a3f8bc6827a55a741a652d34431c7c0fae1a29154dec7710d9bfb4`。
- 作者自检的本次独立运行输出：`/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/position-slab-impl/independent-selftest.json`，SHA同上述bb6b56a0。

完整脚本绝对目录为 `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-r001-audit/planar-centered-square-weak11/computation/`。实际命令见各JSON；只读真实复核接受--module、--input、--slab、--velocities和--output，不调用枚举器大循环。未改main，未提交/推送；一般C₂匹配上界仍未解决。
