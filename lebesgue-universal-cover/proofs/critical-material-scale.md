# 冻结内材料的临界比例

记 $G=120689264254125501/2^{57}$、$a=81827130/2^{28}$、$c=a(1+1/\sqrt3)$。在要求保留第221轮32级和64级的全部母材料、径向材料以及原24个加密原叶的附加锚点时，容许的最小共同缩率是
\[
 s_0=\frac34\!\left(c+\sqrt{c^2+\frac{16}3a^2}\right).
\]
本节从材料居中缩放，所有实际放置的平移和角度仍可分别自由选择。

取 $L=.99999900555134042347742929191955$、$U=.999999005551341$。直接以256位 Arb 检查 $L<s_0<U$。在 $L$ 处，两张表中，除三角材料外的全部八个表与曲体组合在其定义圆盘内均有严格余量。三角部分共有96和192个样点，三角自身的姿态在规范化后固定，核对每个点和三个定义圆盘，只有 $(a,a),(a,-a)$ 对其对应圆盘越界。在缩率 $s$ 处，这一对的包含条件等价于
\[
  2a^2+cs-\tfrac23s^2\le0 .
\]
所以二者恰在 $s_0$ 处取等。各居中凸体含原点，其同心位似关于 $s$ 嵌套。24个加密原叶在 $L$ 处产生的所有内锚点均有严格余量，整数坐标与旧证书使用的逐项相同；径向表也逐项绑定到32级零角行。因此 $s_0$ 对这里列明的完整点表同时必要且充分。[材料根与例外点](../computation/critical-material-scale/records/triangle.json) · [加密锚逐叶记录](../computation/critical-material-scale/records/dense.json)

为进入面积比较，还需排除小面积摆放的原参数域之外、旧纯五体终端以及三十点的支撑帽。取整个 $[L,U]$ 的有向界：对四个曲体中心根域外的最小面积作同径估计，从3,273个原叶重建三角中心投影界，并在全部105条法向区分64条直接材料墙与41条圆盘主导墙。径向裁切给小面积假设下所有材料位于 $\lvert x\rvert<.703$，且双帽仍不交。纯五体终端较高的面积门槛支付缩率损失。原29个帽闭格在整个尺度区间用有理最坏端点核对，覆盖的投影区间也取两端最大。第242轮的有理位移 $p_i-(1-s)z_i$ 的435个点对在全区间保持直径小于 $s$，其2,493行的归一投影损失率与 $s$ 无关。[尺度依赖条件的完整回执](../computation/critical-material-scale/records/interval.json) · [方向投影](../computation/critical-material-scale/records/points/proof.json) · [原面积树第二遍](../rounds/234.md)

对 $s_0<s\le U$，全部所用材料严格内含，可移植原 .83745 全域面积结论。令 $s\downarrow s_0$，
\[
 \Lambda\ \ge\ \frac{G}{s_0^2}
 >0.837451665604544468510297997471 .
\]
同样可改用第244轮缩后的固定十五点：它的原生整数行逐行不低于冻结下端，直径小于 .999448，仍严格小于 $L$。帽公式中实际取同一投影比较，损失为零。此极限也可写成五个连续曲体与一个固定15点输入的下界（最后的单位直径有理组可容纳全部 $s$ 对应的第六体）。[十五点转移](fifteen-point-support-rows.md)

从课题目录调用单入口；TREE 和两个冻结表的来源、散列与整个面积树的重放见第221、234轮。

~~~sh
python -E -B computation/critical-material-scale/verify_interval.py computation "$TREE" "$TABLE32" "$TABLE64" "$WORK/critical" --workers 4
cmp "$WORK/critical/result.json" computation/critical-material-scale/records/result.json
~~~
