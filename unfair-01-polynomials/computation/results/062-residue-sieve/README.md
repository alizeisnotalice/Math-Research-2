# 次数 56 的第 7 项缺口分支精确筛选

[分裂域必要条件](../../../proofs/062-c7-splitting-field-obstructions.md)适用于 \(A_1,B_1>0,\ C_7=0\)，覆盖任意非负实系数、次数和支撑空隙。本次数据为种子 20261118056 所抽到的 163 个次数 56 的多项式。

对素数 \(2,3,7,11,163,1279\) 检验不可约因子次数最小公倍数。151 个目标失败，剩 12 个。Frobenius 幂与多项式 gcd 独立重算 229 项签名。随后用平方自由部分的两处判别式条件重新检查 12 个目标，其模 163、1279 的两处导数 gcd 均为单位。因而这 163 个具体目标全部被排除于上述分支。

部分 Metal 计算在加入第 2 道筛选前已经启动，随后在已运行分段保存后关闭队列，待处理目标不再启动。两类精确结论均来自全整数或有限域运算。

    python3 -B verify_gap7_periods.py inputs-gap7-degree56.json --out /Volumes/data/CodexTemp/T2/a30-gap7-period-replay.json --workers 4
    python3 -B verify_two_ramified_primes.py inputs-gap7-degree56.json --out /Volumes/data/CodexTemp/T2/a30-gap7-ramification-replay.json

运行使用 Python 3.13.18 和 SymPy 1.14.0，生成器使用 NumPy。JSON 含实际版本参数、筛选前缀、目标和已核验回执。
