# 次数 27 和 28 的缺口排除

[完整覆盖和单位理想证书](../../../proofs/064-c7-degrees-27-28.md)。重放覆盖 442,368 和 1,769,472 个早端分支、5,192 和 20,714 个晚端字。主进程对冻结乘子直接重算原式，得到 25,873 个不同的原方程单位恒等式，零余项。

使用 Python 3.13.18、SymPy 1.14.0 和启用断言的解释器。以下命令从本目录运行，输出放在 data 盘临时目录。

    mkdir -p /Volumes/data/CodexTemp/T2/a30-degree28-check
    cp direct_mod_27_28.py modular_27_28.py probe_27_28.py \
       ../060-gap7/verify_degree_19.py ../060-gap7/verify_degrees_17_18.py \
       ../062-gap7/verify_degrees_21_22.py ../063-gap7/verify_degrees_23_24.py \
       ../063-gap7/verify_degrees_25_26.py /Volumes/data/CodexTemp/T2/a30-degree28-check/
    gzip -dc degree-27.json.gz > /Volumes/data/CodexTemp/T2/a30-degree28-check/degree-27.json
    gzip -dc degree-28.json.gz > /Volumes/data/CodexTemp/T2/a30-degree28-check/degree-28.json
    cd /Volumes/data/CodexTemp/T2/a30-degree28-check
    python3 -B direct_mod_27_28.py --only 27 --witness degree-27.json --out degree-27-new.json
    python3 -B direct_mod_27_28.py --only 28 --witness degree-28.json --out degree-28-new.json

原始见证、压缩版本、运行时间和校验哈希见 verification.json，主进程复验摘要见 parent-replay-summary.json。去掉 --witness 会从头生成全部乘子。
