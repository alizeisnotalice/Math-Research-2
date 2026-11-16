# 第34轮解析证明复核

[统一定理](../../../proofs/034-short-factor-isolation.md)在完整交叉替换假设下，关闭所有互素 \(p<q<2p\)、短侧次数 \(<2q\)，以及 \((p,q)=(2,5)\)、短侧次数 \(<10\)。另一侧次数不限。一般均匀和猜想仍开放。

新增一般四角单位引理：若短侧配对支撑为 \(T\cup(T+q)\)、\(t=\max T\)，则任何带符号反设必须有 \(g_0=g_t=g_q=g_{q+t}=1\)。它不确定内部点的系数。

独立终审绑定最终正文SHA，检查倒数模型常数归一、缩放后的正端点、孤立余类、两个固定整数G的全部首负补偿。主线程标准库整数脚本在Python3.13.18实际检查3043个参数余类窗口、320个固定G系数式、1020个倒数支撑变换；另外复算[单滤波抽象反例](../../../proofs/034-abstract-filter-obstruction.md)。这些有限检查不替代解析证明。

    /Users/zzh/.venvs/main/bin/python check_identities.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round034-identity-recheck.json
