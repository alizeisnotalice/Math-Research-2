# A28／A165：二次 Fourier 谱支撑下界

## 问题与假设
均匀测度下 g:F₂^m→{−1,1}，S={α:ĝ(α)≠0}，s=|S|，d=dim span S（线性维数）。目标：证明绝对常数 c>0 使 s≥cd²，或构造 d→∞、s/d²→0 的显式布尔族。原题采用 {0,1} 编码并要求所有一阶系数非零。

## 当前状态
第16轮完成，原题未解决。若旧码C含全部仿射相位并保留各分支独立参数，则满选择器的一切新合法方向仍必须沿用原三角结构；其包含极大扩张恰逐支取旧码的包含极大扩张，最大完成维数为Qm(C)+Q−1、最大完成容量为QΛ(C)。一个指定十一维六变量子码恰有2025个包含极大扩张，全部十五维、容量16；CPU与Metal各核验预测的1605632个合法相位。删减分支后不能沿用此刚性；仅保留独立常数时旧q₀=2归约为恒秩二次型和四点频率分割，联合支撑一般非仿射。耦合共同秩与一般A28仍开放。

研究稿：[第1轮](proofs/001-research-zh.md)、[第2轮](proofs/002-affine-partitions-zh.md)、[第3轮](proofs/003-shared-label-rank-budget-zh.md)、[第4轮](proofs/004-signcode-lifts-zh.md)、[第5轮](proofs/005-matching-and-heavy-normal-form-zh.md)、[第6轮](proofs/006-disjoint-switch-closure-zh.md)、[第7轮](proofs/007-minimum-switch-rank-zh.md)、[第8轮](rounds/008.md)、[第9轮](rounds/009.md)、[第10轮](rounds/010.md)、[第11轮](rounds/011.md)、[第12轮](rounds/012.md)、[第13轮](rounds/013.md)、[第14轮](rounds/014.md)、[第15轮](rounds/015.md)、[第16轮](rounds/016.md)；证据与复现见 rounds/001.md 至 rounds/016.md 和 STATE.md。

## 材料结构
proofs/ 保存证明；computation/ 保存精确计算与最小证书；references/ 保存文献与来源；rounds/ 保存逐轮报告。按 T2 模板组织。

## 工作区
独立工作区 T2-A28，分支 codex/a28-fourier-sparsity；其他课题工作区不修改。临时产物放 /Volumes/data/CodexTemp/T2/a28。
