# H8(x^2)H8(x^3) 最终精确证书 checker：独立复审

被审脚本为
`/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round025-parent/h8_exact/verify_H8_exact.py`，
审阅版本 SHA-256 为
`aa2c976f27bce1a32ce98b26424ebfafa3a0d26c502660d7a02cb0fd4aaf8eb8`。

## 结论

**通过。** 在第 18 轮根幂积平衡引理及本轮 H8 四块/1482 布局覆盖论证作为
明确外部数学输入的条件下，该 checker 足以严格排除所有 1481 个非标准补集
布局；唯一未作负系数排除的 \((E,\mathrm{code})=(6,0)\) 被精确识别为
\(\{H_8(x^2),H_8(x^3)\}\)。常数根集对应的恒定因子另行解析处理。

这是 H8 的这个二元目标的严格有限证书，不是原一般问题的结论。

## 已复核的闭合链

1. **母式/Galois。** 检查器不用缓存结论：有限域因子逐项相乘、不可约性和导数
   最大公因子均以整数重算。型 \([8]\)、\([1,7]\)、\([1,2,5]\) 分别给
   传递、点稳定子中的 7-循环和一个换位；由二重传递及换位共轭得到 \(S_8\)。
   同时固定 H8 的系数、首末项、非自反性和二元 C 的精确卷积。

2. **八个母根和四十个目标根。** 四个正虚部母圆盘及其共轭圆盘均作精确
   Rouché 不等式并两两分离。每个 lift 以其绑定的母圆盘、正确的 \(q=2\) 或
   \(q=3\) 重算 Rouché 上下界；20 个 lift 圆盘及其共轭共 40 个，并全部
   分离。故每个圆盘恰含相应根，40 个根已达 \(C\) 的次数，因而也排除了
   共根/重根遗漏。

3. **标签与根组映射。** 输入群编号必须正好是 0--19；balls 的
   `source_group_index` 亦须一一对应。对每项重验
   `(block,family,lift_index,exponent)`、ID、父圆盘、中心和 Rouché记录；并逐块
   强制两个 a-lift 和三个 b-lift 的槽位、以及输入 block 的相同列表。此前
   “编号调换使 a-mask 选到 b-root” 的伪造路径在这里会触发断言。

4. **独立枚举。** checker 不读 options cache，直接由
   \(3u+2v=E\) 建立 \(E=2,3,4,5\) 的组合选择；E=6 固定首块为全 a-lift
   并枚举其余三块，故补集商计数严格为
   \(81+16+81+1296+8=1482\)。证书记录、唯一标准项与这套独立编码逐个对齐。

5. **圆盘到负系数。** 对每个共轭 lift 根 \(z\)，脚本从有理圆盘重新构造
   \(t^2-2\Re(z)t+|z|^2\) 的闭区间盒。以 Fraction 的端点乘积递推指定
   系数，且每个记录满足
   \[
      \mathrm{stored\;lo}\leq \mathrm{recomputed\;lo}
      \leq\mathrm{recomputed\;hi}\leq\mathrm{stored\;hi}<0.
   \]
   开 Rouché 圆盘被闭系数盒外包，故严格 `hi<0` 对实际根仍有效；不需要 pin
   根积或常数项，正归一化缩放不改变符号。

6. **标准项。** `(6,0)` 逐群等于所有且仅所有 a-lift，补集为所有 b-lift；
   已认证的根数分别为 16、24，配合首一性和相应多项式次数，精确就是
   \(H_8(x^2)\) 和 \(H_8(x^3)\)，而不是浮点近似识别。

## 实际复现

以 `/Users/zzh/.venvs/main/bin/python`、`OMP_NUM_THREADS=1` 运行该脚本，
全程仅整数和 `Fraction` 算术，退出码 0。输出：1482 个规范布局、1481 个严格
负证书、最小认证余量
\(102346591/137438953472\)。独立复跑记录为
`/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round025-parent/h8_exact/H8-a2b3-r8-independent-review.json`，
SHA-256 `500d0bb97d0a22d935e95160d2b1d5424197f8a8898305c6b9f8ad0a27d7720e`。
