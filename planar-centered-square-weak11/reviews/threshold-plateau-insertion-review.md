# 第22轮附加独审：单尺度临界平台

2026-11-15。**审查通过，无数学修改要求。** 本记录独立追加，不改变已锁定的原子跨尺度证明及其独审。

审读附件为作者工作区 `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-cross-scale/planar-centered-square-weak11/proofs/threshold-plateau-insertion.md`，SHA-256 为 `8f313777fde68863c843923da6c5b8cbdb3b859a82f54e44bde5875b7333e6e2`。作者最小脚本及输出的 SHA 分别为 `a817e7ba9f9047b0bd29d6fa9a3686ba6d676192dca1b3d67e3c60d304f7c469`、`5ba35835deb1f58dd308b35617b8c70f5a19ee7d0d7aa29161ff1df875c4902f`；本次按正文独立推导，未把作者的通过标签当成证明。

设 A=[−1,1]²、B=[−1/4,1/4]²，f=1_A、g=4ε1_B、0<ε≤1/4。确有 |A|=4、|B|=1/4、∫g=ε。由于 f≤1、g≤1，**每个中心、每个正半径**的各自平均均≤1，故阈值 1 的严格晕区均为空；ε=1/4 时 g 的等号平台也不属于严格晕区。

对任何 x∈B，取 r=1/2，逐坐标有 B⊂Q(x,r)⊂[−3/4,3/4]²⊂A。这包括 B 的边界和角点：闭见证盒恰好接触 g 的支持边界时仍包含全部支持，积分不丢质量。Q 的面积为 1，所以 ∫_Q f=1、∫_Q g=ε，混合平均恰为 1+ε>1。故整个 B 都是纯混合新增接收区，面积至少 1/4，与新增质量之比至少 1/(4ε)。固定旧输入 f、固定见证半径便已得到无界比值。

上述论证固定严格阈值 1，足以反驳所述统一预算。任意其他阈值 λ>0 可把两个密度同乘 λ 来缩放同一反例；并未声称不变的两个输入在所有阈值下各自晕区都为空。**不能把本例直接改成非严格门槛 ≥1**：此时 B 已属于 f 的非严格晕区，原来那块纯混合新增区域不再被保留。

范围说明准确：不受限 L¹ 预算的失败本身不需要跨尺度机制；原构造的独立内容是有限原子输入加入单个原子时，固定新增质量可以在 N 个分离尺度上产生显式可加的严格混合区域，以及相应定量平滑传递。两者都不推出 C₂ 的新界，也不否定使用旧质量储备的预算。

本次另实际运行一次 Python 3.15.0 `Fraction` 最小核验，exit_code=0；ε∈{1/4,1/1024}，每例检查 B 的 4 个角点、4 个边中点和中心，总计 18 点。直接计算方形交叠长度的乘积，均得 f 积分 1、g 积分 ε；混合平均分别为 5/4 和 1025/1024。没有随机数或浮点容差，未启动大计算或 MATLAB。输出保存在 `computation/results/cross-scale-insertion-audit/threshold-plateau-audit.json`，包含四份最终来源哈希；实际核对原跨尺度 proof 及原 review 的哈希均未改变。

实际运行的核心代码如下，完整有限算术检查可直接复现；全 ε、全接收区域的命题由上面的解析论证承担：

```python
from fractions import Fraction as F
r = F(1, 2)
def length(a, t):
    return max(F(0), min(a+r, t)-max(a-r, -t))
for epsilon in (F(1, 4), F(1, 1024)):
    assert 0 < 4*epsilon <= 1 and 4*epsilon*F(1, 4) == epsilon
    for u in (-F(1, 4), F(0), F(1, 4)):
        for v in (-F(1, 4), F(0), F(1, 4)):
            integral_f = length(u, F(1))*length(v, F(1))
            integral_g = 4*epsilon*length(u, F(1, 4))*length(v, F(1, 4))
            assert integral_f == 1 and integral_g == epsilon
            assert (integral_f+integral_g)/(4*r*r) == 1+epsilon > 1
```

本分工未提交或推送，由主 Agent 合并归档。精确 C₂ 及匹配的一般上界仍开放。
