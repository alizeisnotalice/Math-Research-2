# A28／A165：二次 Fourier 谱支撑下界

## 问题与假设
均匀测度下 g:F₂^m→{−1,1}，S={α:ĝ(α)≠0}，s=|S|，d=dim span S（线性维数）。目标：证明绝对常数 c>0 使 s≥cd²，或构造 d→∞、s/d²→0 的显式布尔族。原题采用 {0,1} 编码并要求所有一阶系数非零。

## 当前状态
第19轮完成，原题未解决。可线性投影的两倍块与最小块组成独立完整分割时，任意共同bent空间满足dimV≤9q/2−1；全公共因子情形为2q−1。部分覆盖仍保留残余维数。完整有限计算覆盖RM(3,6)全部232081500个重量16码字，均可投影，因而覆盖所有六变量真实合法两倍块；不外推高维。MacBook CPU、本机12线程独立算法与MacBook Metal交叉核验通过。一般可投影性、更大块、缺少独立分割的空间及A28仍开放。

研究稿：[第1轮](proofs/001-research-zh.md)、[第2轮](proofs/002-affine-partitions-zh.md)、[第3轮](proofs/003-shared-label-rank-budget-zh.md)、[第4轮](proofs/004-signcode-lifts-zh.md)、[第5轮](proofs/005-matching-and-heavy-normal-form-zh.md)、[第6轮](proofs/006-disjoint-switch-closure-zh.md)、[第7轮](proofs/007-minimum-switch-rank-zh.md)、[第8轮](rounds/008.md)、[第9轮](rounds/009.md)、[第10轮](rounds/010.md)、[第11轮](rounds/011.md)、[第12轮](rounds/012.md)、[第13轮](rounds/013.md)、[第14轮](rounds/014.md)、[第15轮](rounds/015.md)、[第16轮](rounds/016.md)、[第17轮](rounds/017.md)、[第18轮](rounds/018.md)、[第19轮](rounds/019.md)；证据与复现见 rounds/001.md 至 rounds/019.md 和 STATE.md。

## 材料结构
proofs/ 保存证明；computation/ 保存精确计算与最小证书；references/ 保存文献与来源；rounds/ 保存逐轮报告。按 T2 模板组织。

## 工作区
独立工作区 T2-A28，分支 codex/a28-fourier-sparsity；其他课题工作区不修改。临时产物放 /Volumes/data/CodexTemp/T2/a28。
