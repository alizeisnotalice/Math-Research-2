# R57 二次事件处的端点准确比较核

日期：2026-11-15，Asia/Shanghai。本附件只承担二次代数事件处的值、一阶/二阶质量导数比较及分支检查，不承担源目录、bbox 完整性或有效 union 墙结论。R52 冻结文件未修改。

## 1. 有序二次域

固定有理 Q>=0，使用实嵌入 K=Q(sqrt(Q))，sqrt(Q)>=0。元素存成有理对 (a,b)，代表 a+b sqrt(Q)。Q 是有理平方（包括0）时先将 b sqrt(Q) 折入 a。否则 Q>0 且不是有理平方，表示唯一。

对 a+b sqrt(Q) 求符号：a=0 或 b=0 时直接判断；a,b 同号时保留该号；异号时

    sign(a+b sqrt(Q)) = sign(a)*sign(a^2-b^2 Q).

此时才平方，故没有分支伪等号。倒数为 `(a-b sqrt(Q))/(a^2-b^2 Q)`，先拒绝零分母。

在固定事件 t*=alpha+beta sqrt(Q) 处，原始端点

    e(t)=c+s sqrt(q0+q1*(t-t0)), s in {-1,0,1}

化成 `c+s sqrt(U)`，其中 U=q0+q1*(t*-t0) 属于 K，U>=0。c 可在底层比较核中属于 K；路径 recipe 接口的 c、q0、q1 均为有理数。位置速度单独保留在 recipe 元数据中。

## 2. 一根式与两根式比较

令 h,k,U 属于 K，U>=0。准确符号函数 S(h,k,U)=sign(h+k sqrt(U)) 采用以下分支：

- U=0 或 k=0：返回 sign(h)。h=0 且 k!=0、U>0：返回 sign(k)。
- h,k 同号：返回该号。
- 仅当 h,k 异号时：返回 `sign(h)*sign(h^2-k^2 U)`。

最后一支已确认两项相反，平方后的零正好对应真实等号。

比较 e1=c1+s sqrt(U)、e2=c2+r sqrt(V)，令 d=c1-c2、y=d+s sqrt(U)，先用 S 得到 sign(y)。若 r=0 或 V=0，直接返回 sign(y)。若 y=0，则返回 -r；若 y 和 -r 同号，则返回 sign(y)。其余情形 y 与 -r sqrt(V) 异号，故

    sign(e1-e2) = sign(y)*sign(y^2-V),
    y^2-V = d^2+s^2 U-V + 2ds sqrt(U).

再调用一次 S 即完成准确比较。整个过程只在 K 内加、乘、除和判断符号，且每次平方前已有对应的符号条件。没有浮点容差，也不以区间重叠推出等号。

## 3. 持久身份、非持久等号与左右序

原始 recipe 是 `[c,w,s,q0,q1]`；w 为物理 z 位移速度，不参与本节 t 导数。持久身份须保留整个根式轨迹，而不能仅比较 t* 处的 q 值。有理平方只在 q1=0 时折入常数；q1!=0 时即使 q0 为有理平方也保留。不同 t0 的 recipe 用全局仿射常数 q0-q1*t0 对齐。

对于 U>0，质量导数和二阶导数分别为

    e' = s q1/(2 sqrt(U)) = [s q1/(2U)] sqrt(U),
    e'' = -s q1^2/(4U sqrt(U)) = [-s q1^2/(4U^2)] sqrt(U).

任意 C sqrt(U) 可准确转成 `sign(C)*sqrt(C^2 U)`，因此仍交给同一端点比较核。q1=0 或 cap 常数支的两阶导数为0；非恒定根式在 U=0 时拒绝解析 jet 认证。

两个有理 recipe 在 t* 的值相等后，先比较一阶 jet；非零差决定右序，左序反号。若一阶相等，比较二阶 jet；非零差决定两侧同序。若两阶都相等，则它们沿 t 持久相等。证明如下：若共同一阶导数为0，两者均为常数；若非零，由 `e''/e'=-q1/(2U)` 相等，得到两者 q1/U 相同。再由一阶导数相等，得到符号相同、U 相同，进而 q1 与常数偏移相同。因此非持久解析接触阶最多为2。

沿 t 持久相等不代表完整二参数身份相同：两个端点可能在固定 z0 具有同一 t 轨迹而 w 不同。输出分别记录 `persistent_in_t` 与 `same_full_recipe`，不把二者混为一谈。

## 4. 已实际执行的证据

标准库 Fraction/isqrt，CPython3.15.0、macOS27.2 arm64，无随机种子，无浮点容差。测试脚本 exit0，完成14个分支边例、150个底域符号/倒数例、900个可还原为 K 元素的平方根比较，以及180个通用嵌套根式比较。通用非零比较用另一套有理定向根区间独立分离；区间不能分离时不会假定相等。另准确核对一阶穿越、二阶相切、质量持久但位置速度不同、参数原点变更，以及四类非法输入拒绝。

实际 R52 两盒墙的最小输入也已独立执行：Q、alpha、beta 见 `inputs/round057-quadratic-kernel/wall-pair.json`。准确检查 qB(t*)=(D-sqrt(Q))^2，并确认 D-sqrt(Q)>0、0<t*<1/2、t* 严格落在冻结128位证书括区内。A.left 与 B.right 在事件处准确相等；A.left-B.right 的一阶导数为正，故左侧 A.left<B.right，右侧 A.left>B.right。它仍只是原始端点穿越，不证明实际并边界暴露。

## 5. 几何批任务接口

模块 `r057_quadratic_endpoint_kernel.py` 提供 `Field`、`Endpoint`、`compare_endpoints`、`recipe_at_event`、`prepare_recipe`、`compare_prepared`、`contact`、`order_cluster`。

CLI 输入为

    {"Q":"有理数", "event":["alpha","beta"], "t0":"1/2",
     "endpoints":[{"id":"标签", "expr":["c","w",s,"q0","q1"]}],
     "candidate_clusters":[["标签1","标签2"]]}

调用 `python3 -B r057_quadratic_endpoint_kernel.py input.json --output out.json`；`--cluster-index N` 可将已生成的独立 x/y 端点簇分配给操作协调者。每簇采用比较排序和相邻比较生成中心等价类及左右序，没有构造全对比较图。输出明确只认证给定簇内部；簇的覆盖完整性和簇间分离由几何生产者另证。

精确源文件 SHA、大小、输入及输出以配套 integration-manifest 为准。主 Agent 负责独立审读和最终整合；本助手未提交或推送。
