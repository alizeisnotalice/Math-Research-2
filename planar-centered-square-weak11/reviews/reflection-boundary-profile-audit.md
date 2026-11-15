# 第37轮多接缝条件引理与平台障碍：Astra独立审读

日期：2026-11-15。审读者 astra_reflection。本稿从两份新证明的公式逐项复核，绑定下列最终版本；旧第30轮文件未改。有限算术输出已实际读取，但不把Sol运行记作本审者重跑。

## 多接缝条件引理

通过。每个候选的真实反射混合横半宽是kr_j−g_j，纵区间是旧区间向两端延伸kr_j。E_j必须避开由全部旧接收矩形构成的V(w_j)，包括右端等号，不能只排除候选自身。

同轴patch横截面为嵌套区间，并长度准确为2max w_j；因此一维端点扫描给真实并面积，不是近似网格。纵重数q控制接收重叠，源占用重数是另一对象。条件α、β、q成立后，下界对所选源质量之和合法，进一步用源取并质量只会减小右端，不需要源彼此不交。

有限L¹传递把全部旧矩形和全部patch用同一个η内缩；并集面积损失≤各矩形损失之和，即使patch原来重叠、部分内缩坍缩仍成立。保留严格阈值的中间半径存在。分离行例D=(96/5)N和η=1/120、增益1/25与平滑后1/50均重新核对。

本引理没有把所有靠边界的源质量自动转成有效patch质量；其外露长度α和接收重数q必须另外认证。将有效质量按好patch定义之后再断言它占全部活跃质量固定比例，将是未证步骤。作者已明确保留此缺口。

## 密度2平台与完整halo上界

通过。矩形外或边界中心在至少一个坐标方向最多捕获方形的一半，因此密度2时全半径平均≤1；内部有小半径严格平均2。这准确给旧halo与反射拼接后的完整halo，而非局部見证估计。

网格量化对所有半径满足μ_h(Q_r)≤2|R∩Q_(r+h/2)|。中心距R某一坐标至少h时，右端≤4(r+h/2)(r−h/2)_+≤4r²，严格排除全部外侧中心。反射后两份网格的参考方格仍不交、密度至多2，全部位于B_h，故同一论证给完整新halo外界。真实新增面积≤4hL+8ah+8h²是上界，不是由某个小patch的下界推断。

## 准确最右端owner及高比值核心

通过。任意bbox跨p列q行，质量≤2h²pq且正纵宽强迫(q−1)²<2pq，整数性推出q≤2p+1。右端上界的根式差√(p²+p/2)−p严格随p增加。因此全局最大必须在p=2n、q=4n+1的满bbox取得，且质量上界必须取等。这完整覆盖非饱和任意原子子集。

K≥2保证这些滑动行块存在。每个块纵跨度4a，临界半径√(4a²+ah)在2a与2a+h之间，选近临界但更小的半径能完整捕获该块、排除邻行。全部准确最右owner的源并为整个网格，最右且最高owner也有不退化质量16a²+4ah。由真实halo外界，它的质量仍不能单独保证统一正增益。

核心隔离使用全输入总质量T和源间距>√(2T)，所以反射前后任何严格见证均不能同时捕获不同组件的正质量；对L¹部分捕获和有限原子都成立。两份核心与平台的halo亦分离，真实接缝和新增面积因此由平台控制。

高比值量词保留为：对每个固定q<C₂，可选一个核心及固定正平台质量m，使J>q，且沿网格细化m/T保持正值而增益趋零。m/T可依赖q；本附件没有声称q趋C₂时仍有统一正下界，也没有给出同一固定有限原子数的近极值列。

## 实际结果和结论范围

两份Sol输出均为Python3.15.0 Fraction、单进程、exit_code=0的实际结果。patch例覆盖N=1,2,5,20,64及重叠/遮挡检查；平台例共7个(K,n)，最大512原子，仅检查准确bbox形状、质量覆盖、全半径标量界，没有声称计算完整512原子halo。一般结论依赖上面的解析论证。

这两份附件与竖格完整面积正例共同说明：外露边界轮廓可以产生线性累计，但准确边界owner的质量不足以保证该外露量。当前高比值周期模板的边界改善仍未认证；C₂精确值未解决。正式提交推送由root负责。

## 版本绑定

- `pcw-sol-reflection-geometry/proofs/reflection-multiple-patch-packing.md`：10419 bytes，SHA-256 `a03dd17c8efee63a7ae204f604aca3b32a8b2932f3ec9422932693fb51b6c334`。
- `pcw-sol-reflection-geometry/computation/reflection_patch_packing_check.py`：3248 bytes，SHA-256 `0745587881b6854f450a0df1906ffe552e7f9ff195c190b4a907b7dbb1c028fb`。
- `pcw-sol-reflection-geometry/computation/results/reflection-patch-packing-check.json`：1919 bytes，SHA-256 `3c51e4804071f1ef39571638bd27470c23c24de803fc5b9dd50dadd09d32215d`。
- `pcw-sol-reflection-limits/proofs/seam-budget-plateau-obstruction.md`：14591 bytes，SHA-256 `789e8172c5659a693a821bc47125182e6fb770287996efc36c3f932fba8df40b`。
- `pcw-sol-reflection-limits/computation/seam_budget_plateau_obstruction.py`：5329 bytes，SHA-256 `64ec77f7d769c8d0505367ee92f6d60a6db15b7fbe71799be7700e150e15a6f4`。
- `pcw-sol-reflection-limits/computation/results/seam-budget-plateau/exact-check.json`：5279 bytes，SHA-256 `97f340848b6b4758210e9fba06671ddebfa4e57a03e72985eafdfb6f95300870`。
