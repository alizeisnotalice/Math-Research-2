# 第 7 项缺口与次数 20

[中文证明](../../../proofs/061-c7-degree-20.md)。假设归一化的两因子一次项均严格正，且 \(C_7=0\)。任一因子次数为 20 会使第 20 项的内部卷积项全部为零。两侧第 8、9、10 项的 18 种占用方式中，10 种给非整范数，4 种给卷积预算超支，4 种由模 163 的 512 个二元字排除。

使用 Python 3.13.18 标准库和启用的断言。将三份脚本放到同一个 data 盘临时目录，再运行：

    mkdir -p /Volumes/data/CodexTemp/T2/a30-degree20-check
    cp verify_degree_20.py ../060-gap7/verify_degree_19.py ../060-gap7/verify_degrees_17_18.py /Volumes/data/CodexTemp/T2/a30-degree20-check/
    python3 -B /Volumes/data/CodexTemp/T2/a30-degree20-check/verify_degree_20.py --out /Volumes/data/CodexTemp/T2/a30-degree20-check/new.json

主进程的实际回执为 [degree-20.json](degree-20.json)。重放与首次验证的数学字段一致。
