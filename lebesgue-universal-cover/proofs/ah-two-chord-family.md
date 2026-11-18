# 四个参数共用一个摆放的双线宿主

在直线 Ta–Z 上取横坐标 w，向宿主内取纵坐标 y。铰点 p、转角 δ 的保留半平面为 y≥tan(δ)(p−w)。两条直线和原边界给出的保留高度是三项 0、tan(δ₁)(p₁−w)、tan(δ₂)(p₂−w) 的最大值。参数独立地取遍

| 裁线 | 铰点 | 转角（弧度） |
|---|---|---|
| 尾线 | [0.4268,0.427] | [0.00378,0.0038] |
| 前线 | [0.394,0.395] | [0.00480,0.00483] |

对任何实数 w，截断高度 max(0,tan(δ)(p−w)) 对 p、δ 分别单调不减。因此四个参数同时取上端时得到最小宿主，同时取下端时得到最大宿主。最小宿主中对给定输入认证的一个摆放同时适用于全盒。最小和最大宿主的真实边界分别按 Green 面积公式用有向区间积分。

根分支还可转移严格余量。记 b=tan(δ)，以加号标出参数盒上端，则
\[
y+b(w-p)=\frac b{b_+}\{y+b_+(w-p_+)\}+(1-b/b_+)y+b(p_+-p).
\]
末两项非负，故有符号距离下界至少保留最深角点的 b/b₊ 倍，法向归一化的附加因子不小于一。两个切面分别使用此恒等式，取两者较小余量。真配对圆弧和负 w 段的几何估计、两次实际命中按亏损选首点、同一根处的其余边界判定，在最深角点一并重放。第191轮的十一段 EH 裁面权重在该角点重新校验。扣除共同代价后，整族跨宿主摆放的净余量仍大于 0.000001079977。

以下命令在课题目录运行，三个输出可以与同目录的 first-hit-family-certificate、first-hit-family-eh、first-hit-family-combined 三份 JSON 逐项比对。

```sh
python -B computation/certify_ah_dual_family.py computation computation/first-hit-family-repaired-input.json AH.json
python -B computation/certify_eh_upper_chain.py computation AH.json computation/eh-upper-chain-parameters.json computation/eh-upper-chain-weights.json EH.json
python -B computation/combine_ah_family_eh.py AH.json EH.json family.json
```
