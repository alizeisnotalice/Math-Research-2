# 0.006 单弦同十一段 EH 的直接支撑复算

保持 Gibbs 角 (1.5497^circ)。AH 线取铰点 (p=0.427)、角 (0.006) 弧度。AH 单弦作用域、定宽双命中几何和 Green 帽面积沿用同一份冻结单弦证书。EH 使用已经验证的十一段轮廓、十七种有序命中条件和63份非负对偶。

组合时先从原参数重算最深包络 EH 弦与旧包络弦的局部替换代价。各 AH 参考格余量有足够预算时，直接减去此代价。其余格属于共同角括编号 0、1、30。重新构造原 AH 图，在派生几何 SHA-256 完全对应的前提下，加入由最深 EH 上端和旧下端确定的精确单位支撑行，再验证这三个闭格的相容根和同一姿态。格0、1沿用按同一输入支撑值作出的两根分支，格30沿用原高根。所有验证同时要求原有边及新增正余量。十一段实际耳落在包络单弦所裁耳内，因此这个新增真支撑行充分。AH 与 EH 的裁切面积可相加，因为 EH 耳的整个外包矩形严格通过 AH 裁线。

冻结单弦证明与十一段权重合并归档后，可从课题目录重放：

```sh
python -B computation/certify_eh_on_flexible_chord.py computation computation \
  computation/ah-flexible-chord-certificate.json \
  computation/eh-upper-chain-parameters.json computation/eh-upper-chain-weights.json JOINT.json
```

第一次复算把两个 `computation` 参数分别指向主工作区的 EH 和独立 AH 工作区。程序对两份原输入、派生图及最后 EH 变体均记录指纹。
