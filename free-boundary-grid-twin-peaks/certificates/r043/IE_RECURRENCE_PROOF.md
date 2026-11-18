# 单区域根带一峰计数：独立集容斥递推

固定 (P_n\square P_n) 的诱导区域 (H)，并令 (A\subset H) 是允许根集：
根在外根带内，且它在全图中的邻点均在 (H)。标号按从大到小揭示。对
一个独立集 (S\subset H)，令 (N_H(S)) 为其在 (H) 中的开邻域，并令

\[
C_H(\varnothing)=|H|!,\qquad
C_H(S)=\frac{\sum_{v\in S}C_H(S\setminus\{v\})}
 {|S\cup N_H(S)|}\quad(S\ne\varnothing).
\]

这里 (C_H(S)) 是全部 (|H|!) 个标号中 (S) 的每个顶点都是局部峰的
精确个数。令 (T=S\cup N_H(S))。当 (S) 独立时，(T) 中最大标签必在
(S)。条件于该最大标签为 (v\in S)，其概率为 (1/|T|)，且其余顶点的
相对顺序仍均匀。又因 (S) 独立，(v) 不在 (S\setminus\{v\}) 的闭邻域，
故剩余局部峰约束的概率正是 (C_H(S\setminus\{v\})/|H|!)。对所有 (v)
求和并乘 (|H|!) 即得递推。因此该除法逐步整除；程序对每一步余数作
`assert`，不使用浮点数。

令 (B_v) 为「(v) 是 (H) 的局部峰」；共同局部峰只能形成独立集。于是
逐个实际峰集的容斥给出允许根带中的唯一局部峰数

\[
a_H(A)=\sum_{\varnothing\ne S\text{ independent}}
 (-1)^{|S|-1}|S\cap A|C_H(S).
\]

对轮廓 (h=(h_y))，取 (H=L_h=\{(x,y):x<h_y\}\)。右区域经左右反射与
左轮廓 (n-h) 同构，因此其根带一峰数为 (a_{L_{n-h}})，但程序实际分别
计算 (h) 与 (n-h)，不把面积不同的两者误认为相等。两侧标签交错的精确
事件数为

\[
C_h=\binom{n^2}{|L_h|}\,a_{L_h}\,a_{L_{n-h}}.
\]

本轮只计算 (n=7\), (h_y\in\{3,4\}) 的 128 条曲线。每边至多 28 个顶点、
独立数至多 14；正、负容斥和的粗界
(28!\cdot14\cdot2^{13}<2^{128})，故区域递推与容斥使用无符号 128 位整数。
最终 (C_h) 另以 256 位整数输出，避免把两个区域计数及 shuffle 强行压入
128 位。这个精确有限计算不控制 (n\to\infty) 的归一化质量。

## n=5 门槛与 n=7 归一化来源

门槛逐一比较全部 (h_y\in\{2,3\}) 的 32 条曲线，与
`/Volumes/data/CodexTemp/T2/a022-r042-diagonal.stdout` 的
`left_rooted_count`、`right_rooted_count` 和 `count` 相等。

若需要 (n=7) 的全图双峰归一化，本轮只引用公开项目的精确整数记录
`/Volumes/data/CodexProjects/T2/free-boundary-grid-twin-peaks/computation/results/n7-exact-r005.json`
(SHA-256 `ff7974ff8487a339516fee95324674efd873abd9a36b8fc50b0f8121c632f8d6`)，
其 `A2` 字段为
`137077928808481869759712347602605837601908319657172480`。这不是 SMC 输出。
