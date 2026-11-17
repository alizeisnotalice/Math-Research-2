# A28／A165：二次 Fourier 谱支撑下界

## 问题与假设
均匀测度下 g:F₂^m→{−1,1}，S={α:ĝ(α)≠0}，s=|S|，d=dim span S（线性维数）。目标：证明绝对常数 c>0 使 s≥cd²，或构造 d→∞、s/d²→0 的显式布尔族。原题采用 {0,1} 编码并要求所有一阶系数非零。

## 当前状态
第24轮完成，原题未解决。一般共同bent空间可在O(k²)余维内删去全部<2q方向；显式次数截断三角子码满足D=2q−1−2Σ_{d<ℓ}C(k,d)、最小重量2^ℓq。可使D/q→2同时最小重量/q→∞，因此无条件低重量覆盖路线不足；这些子码并非包含极大，不是A28反例。

研究稿：[第1轮](proofs/001-research-zh.md)、[第2轮](proofs/002-affine-partitions-zh.md)、[第3轮](proofs/003-shared-label-rank-budget-zh.md)、[第4轮](proofs/004-signcode-lifts-zh.md)、[第5轮](proofs/005-matching-and-heavy-normal-form-zh.md)、[第6轮](proofs/006-disjoint-switch-closure-zh.md)、[第7轮](proofs/007-minimum-switch-rank-zh.md)、[第8轮](rounds/008.md)、[第9轮](rounds/009.md)、[第10轮](rounds/010.md)、[第11轮](rounds/011.md)、[第12轮](rounds/012.md)、[第13轮](rounds/013.md)、[第14轮](rounds/014.md)、[第15轮](rounds/015.md)、[第16轮](rounds/016.md)、[第17轮](rounds/017.md)、[第18轮](rounds/018.md)、[第19轮](rounds/019.md)、[第20轮](rounds/020.md)、[第21轮](rounds/021.md)、[第22轮](rounds/022.md)、[第23轮](rounds/023.md)、[第24轮](rounds/024.md)；证据与复现见 rounds/001.md 至 rounds/024.md 和 STATE.md。

## 材料结构
proofs/ 保存证明；computation/ 保存精确计算与最小证书；references/ 保存文献与来源；rounds/ 保存逐轮报告。按 T2 模板组织。

## 工作区
独立工作区 T2-A28，分支 codex/a28-fourier-sparsity；其他课题工作区不修改。临时产物放 /Volumes/data/CodexTemp/T2/a28。
