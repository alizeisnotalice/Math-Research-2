# R59 新边界质量手术算核独立审查

2026-11-15；sol_reflection_geometry。审读作者 `computation/boundary_weight_surgery.py`，实际字节SHA256 `1a0233317706618b5f2d1231d10ac7da1e5f46324e1bacc9320eb4ac9632043d`。结论：t=1/2 两侧实际轴、source/setup、完整面积 loop 和真轴误差接口通过；可执行两侧完整面积。本审查没有重跑任何大型halo面积。统一于t∈[0,1]、N>=7的解析合法性见 `proofs/r059-surgery-uniform-tail-locality-independent.md`，不由有限样点冒充全族证明。

## 1. setup/source和退化端点

left使用不取模x'=X-x，I0列0 factor=t，其他标签列0 factor=1；right保持alpha，I0列0 factor=1-t，其他标签列0 factor=0。负列全部full。恰对应定义的两局部半平面。actual xmax逐标签取z>0时a、z=0时a-Px，再取全标签max，故不误把已删除的col0当正质量source。

新独立脚本 `computation/r059_surgery_axis_kernel_independent.py` 自行写setup与inclusive格点枚举，实际对t=0,1/2,1两侧全部source字典（含精确质量）与被审函数逐项相等；zero列端点覆盖。target导入时禁用bytecode，只读别人工作区，不调用大型area。有限样点的字典比对只是实现测试，全参数正确性来自setup及列权重公式的直接核对。

## 2. 完整轴目录和独立根式

轴bootstrap使用实际xmax，逐标签枚举所有 a+iPx>xmax-R、i<=0 的正质量source，nr=ceil(2R/Py)+1，B为全部row质量乘nr；rootup给4Rnext²>B。每轮row标签/质量与独立重建准确一致。任一有竞争力的右端u>xmax，其owner最左源>xmax-R，因此这个目录完整；不是只检查猜测singleton。

作者两轴各48,323个预剪candidate bbox，最终R=543/250，156原子、22x24坐标。本人不使用其dyadic-pruning算法，对两目录各全部75,900坐标bbox独立比较：left5796非空、5785严格、11平手；right6175非空、6164严格、11平手。所有比较都是Fraction和单根式平方判据，正L平方分支次数0。

令q=319996899999/2560000000000，准确得到变换左轴X+sqrt(q)、原左轴ell_(1/2)=-sqrt(q)，原相对NPx右轴r_(1/2)=sqrt(q)。两个作者dyadic区间均由准确平方验证包含该真根。证据保存 `results/round059-surgery-geometry/axis-kernel-review.json`。本机Python3.15.0，实际exit0，无随机数/浮点容差，没有N扫描。

固定64bit脚本在本次t=1/2通过；不能据此声明对任意极小正端点质量都自动成功。若部分列singleton半径小于dyadic分辨率，`lo>xmax`可能安全断言失败，需要增大bits或更强下界；这不会产生错误证书。数学全族尾界与局部性不受该实现精度问题影响。

## 3. 完整面积循环

area的接收窗[anchor,C]x[0,Py]，left anchor=X-3Px，right anchor=-3Px；准确对应原坐标的两个完整半halo窗。source_atoms对原侧用[anchor-R,C+R]，镜像预像[C-R,2C-anchor+R]，标签逐列质量m*(z if i==0 else1)，正确。

该循环与冻结R46完整bbox逻辑一致：全部source坐标prefix质量，span²>=M的空接收准确排除，cap15后裁切，整数isqrt两向根界，inner为空不丢outer，最后矩形并面积不重复加。scale包括anchor、C、Py和全部坐标分母；R=15整数、lambda=1已断言，质量公分母W独立。所用统一单缝R15尾界已在新理论稿证明，因此是全尺度面积，不是截断算子。此为代码审读，并未执行area函数。

## 4. fixed-anchor真轴扰动：必须先居中

物理坐标中mirror移动2delta、窗口上端移动delta，不能简单声称所有物理端点仅移动delta。正确证明：将每个轴c的坐标都移到x-c，原source移动-delta、mirror移动+delta，接收窗变为[anchor-c,0]，其左端同样移动-delta。lo>实际xmax保证所有source在轴两侧的全排序固定。每个固定索引bbox的质量和cap-root不变，裁切后的左右接收端点各为若干1-Lipschitz函数的max/min，仍各至多移动delta。因此每bbox面积对称差<=2delta Py；包含激活、空盒、移动裁切墙。并集差由全部bbox之和支付，面积在轴居中平移下不变。

共同目录确实可取低轴lo：原source下cut anchor-R固定，上cut c+R高于全部原source不生效；mirror下cut c-R随c增大而收紧，上cut2c-anchor+R高于全部原source不生效。独立比较高轴标签集是低轴标签子集，并逐项比对作者source字典。

两侧共同目录均9984原子、423x85坐标，K=327765780全部bbox；每侧准确误差

2delta Py K=2877840232728099021017613/14757395258967641292800000000000000。

这包括moving anchor-c，不应额外再补一项Py*delta，也不能将当前active盒数代替K。实际证据输出保存两个catalog和轴SHA。后续最终area输出须绑定本kernel字节与本轴输入；如kernel再改须保留版本映射，不能覆盖历史来源。

## 5. 尚待面积

本次仅确认可合法计算Lambda_L(1/2)、Lambda_R(1/2)，然后严格组合K=Lambda_L+Lambda_R-6A。两根恰对称不意味着两个全halo面积相等；原template左右内部标签仍不同。K(1/2)符号、新下界、其他t最优性与一般C2均未判定。本人无提交推送；root统一整合。
