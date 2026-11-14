# 第4轮分组位置与质量联合调整：独立完整认证

冻结输入 `group-search/epoch7/candidate_0504.json`，SHA-256 `e1f4cdb0a1868645677f306777c6f1368d41919c63d56fde530df809f1de8c42`。沿用原独立面积程序（SHA-256 `d1e05a0627d699b54472b1b3e74ad500035135853298d33f519f78d3a17903e8`），重新生成全部窗口原子、坐标框和内外并集，没有读取主计算实现。

## 输入变化与合法性

与第3轮最终准确输入逐字段使用 `Fraction` 比较：P_x、P_y 和 λ=1 均不变；位置和质量都有调整，并非纯质量扰动。共有4个 x 坐标、24个 y 坐标及27个质量变化。51个新质量全部严格为正，零/负质量均为0，最小质量为 `1/1250`；全部新位置仍在基本半开胞元内。

准确总质量仍为 `65941670400277227073796193886781311/12500000000000000000000000000000000`，总质量差为0。27个质量变化中1个增加、26个减少；总转移质量为 `9061572822506141189133794218655341/375000000000000000000000000000000000000`（约 `2.41641941933e-05`），质量 L¹ 差为 `9061572822506141189133794218655341/187500000000000000000000000000000000000`。全部准确变化另保存在证书中，分组核对如下：

| 原子/字段 | 相对第3轮的准确变换 |
| --- | --- |
| 23、24 的 x | 减 P_x/100000 |
| 25、26 的 x | 加 P_x/100000 |
| 27、29、…、49 的 y | 减 P_y/10000 |
| 28、30、…、50 的 y | 加 11P_y/100000 |
| 23、24 的质量 | 乘 999/1000 |
| 27–50 的质量 | 乘 9999/10000 |
| 2 的质量 | 增加全部被移除的质量，准确增量 `9061572822506141189133794218655341/375000000000000000000000000000000000000` |
| 其他字段 | 不变 |

## 尾界、q与完整几何

从新输入准确重算 `g(5)=m/(P_xP_y)(1+P_x/10)(1+P_y/10)`，结果为

`23820787253964867284592995960140706111190198205884874124835164741/24052915052186143128797603913750000000000000000000000000000000000`。

严格小于1，余量为 `232127798221275844204607953609293888809801794115125875164835259/24052915052186143128797603913750000000000000000000000000000000000`；因总质量和周期确实未变，重算值与第3轮相同。没有仅凭历史尾界继承全半径结论。

q准确重算为 `3`。窗口实际 k、l 指标范围分别为 `[-2, 2]`、`[-2, 2]`，全部包含于[-q,q]；每原子指标范围重新计数给出1016个副本。

以 bits=40、cap=5 完整枚举2,329,036个跨度合格坐标框，其中28,019个零质量框（不含原子的几何框，并非零质量源），2,301,017个正质量框；内外各6,810个裁剪后正面积候选，各6,397个去重矩形。18个小例及严格根号门槛、cap、裁剪、外正内空测试重新通过。面积以独立逐条带活动区间合并精确计算。

## 准确包围与严格改善

内面积为 `482260334994160682089760129500500502594451975530031442985273/53687091200000000000000000000000000000000000000000000000000`，外面积为 `1446781004988841165288195492155348654566618036346564055518319/161061273600000000000000000000000000000000000000000000000000`。

周期弱比值下界为 `482260334994160682089760129500500502594451975530031442985273/283217317812801919615520431313000857449005056000000000000000`，上界为 `1446781004988841165288195492155348654566618036346564055518319/849651953438405758846561293939002572347015168000000000000000`，阅读用小数约为 `[1.7027925365528678, 1.7027925365603522]`。严格性由准确分数承担。

| 与第3轮的差 | 准确正值 | 阅读用近似 |
| --- | --- | --- |
| 第4轮下界 − 第3轮下界 | `11959238729193398647441484289836342125484360183851/5437772502005796856617992281209616463020897075200000000` | `2.19929000795e-06` |
| 第4轮下界 − 第3轮完整外界 | `47836782807927354583709239646201103926414976913529/21751090008023187426471969124838465852083588300800000000` | `2.1992820953e-06` |

两个正号均用准确 `Fraction` 差分断言，所以严格改善超过完整认证误差。位置和质量同时改变，本次对冻结的新准确输入已经重做完整几何与面积，没有以旧几何或差分导数代替。

已审有限 L¹ 传递仍适用：q=3、D=2^40、ε=2^-41，每个整数 N>6 对应输入的有限比值至少 `(1−6/N)²` 乘本次下界。一般最优性、匹配上界和该联合族的全局优化均未证明。

## 实际运行与复现

实际运行 UTC `2026-11-14T16:51:32.516544+00:00`，Python `3.9.6 (default, Oct 19 2026, 19:18:13)`，平台 `macOS-27.2-arm64-arm-64bit`，完整枚举与面积耗时 `4.147` 秒。结果 `/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round004/independent_final.json` 保存全部准确结果、输入变化、tail/q复核、正差和版本摘要。

实际命令（独立工作区根目录）：

```sh
TMPDIR=/Volumes/data/CodexTemp/T2 TMP=/Volumes/data/CodexTemp/T2 TEMP=/Volumes/data/CodexTemp/T2 XDG_CACHE_HOME=/Volumes/data/CodexTemp/cache PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 planar-centered-square-weak11/computation/independent_periodic_halo.py --input /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round004/group-search/epoch7/candidate_0504.json --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round004/independent_final.json --bits 40
```

没有改动主仓、其他工作区或独立脚本，没有创建提交与推送；由主Agent整合归档。


## Independent reread of the final round-4 proof and report

Actually read `proofs/round004-lower.md`, SHA-256 `43a271601c4bacde205e9761b74df306a59f18f5765a486b10b0b4b9dc427b48`, and `rounds/004.md`, SHA-256 `868465eb3bb15be624277abeea6056dbcd24cb3e866e0a1b4d1bb51ac1747297`, at UTC `2026-11-14T17:05:15.259861+00:00`. No geometry was rerun.

PASS: both exact periodic-ratio endpoints, the exact round-4 minus round-3 lower-bound difference, and the lower-bound fraction in the round report agree with the independent certificate. The archived exact input is byte-identical to the frozen candidate. The reported changes (4 x coordinates, 24 y coordinates, 27 masses), positive masses, conserved total mass, and unchanged periods are consistent with the recorded field-by-field check. The round-4 L4 notation explicitly refers to T2 and is distinguished from the legacy L4.

PASS: the finite L1 transfer proves C2 >= L4. The newly recomputed tail is below 1, q=3, and the reviewed round-1 transfer proof is unchanged. The positive finite source, cell coordinates, D=2^40 inner-radius margin and common epsilon=2^-41 meet its hypotheses. For every integer N>6, the smoothed compact input has mass N^2 m and valid receiver area at least (N-6)^2 A-minus. Fixing this exact input and D before N tends to infinity yields the stated non-strict lower bound. The input's periodic upper endpoint is not a general upper bound.

The pronoun correction is present: the subject exceeding the round-3 full outer bound is explicitly the new lower bound L4. No wording issue remains.

This textual audit does not reverify the 4104 search logs, CPU statistics, remote commit identity, or inherited general upper-bound kernel. Those retain their separate execution/provenance scope.
