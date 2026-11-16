# 短G必要支撑模式：完整有限枚举

在带符号反设下，已证 \(\operatorname{supp}G=T\cup(T+q)\)，其中 \(0\in T\subseteq[0,q-1]\)。对这个支撑，\((a,b)\mapsto pa+qb\) 单射等价于不存在差p：若发生碰撞，互素性给 \(a-a'=qk,b'-b=pk\)；支撑宽度小于2q，非零k只能是±1。反向用支撑中的0、q构造碰撞。该判据也是T不含差p或q−p。

若G有模p孤立余类，完整交叉A直接隔离全部F系数。因此本次枚举仅保留没有孤立余类的必要模式，不为它们构造系数。

全范围为互素 \(1<p<q\le20\)，108个参数对。兼容二元支撑加和的T共106708个；18734个模式通过无孤立余类筛选，分布在45个参数对。63个q<2p参数对均无幸存，与解析定理相符。这里的“幸存”只指局部必要规则；其中(2,5)两个模式也已由本轮系数论证全部关闭，不是反例或未认证分解。

生产程序按T的差禁配DFS。独立验证器不导入生产程序，按所有奇数q位掩码枚举，重建每个具体T和G支撑，并检查所有参数对、计数、无重复和无遗漏。Agent原版实际2.652秒；主线程给验证器添加输入/输出参数后，从压缩包解压的输入实际全量重跑3.281秒。两者规范流SHA一致。删除一个幸存模式并同步修改顶层/行计数的生产篡改被拒绝。

原JSON约6.69MB以确定性gzip归档，原始和压缩哈希见compression-manifest.json。解压到data临时目录并核对原始哈希后可运行：

    python3 audit_short_g_singletons.py --input /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round034-parent/short-G-singletons-q20.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/short-G-recheck.json

这些有限支持模式不决定任意实系数可行性，不扩大原均匀和猜想的严格有限次数范围。
