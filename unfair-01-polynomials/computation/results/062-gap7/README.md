# 第 7 项缺口的次数下界与尾部消元

- [次数 21、22 排除](../../../proofs/062-c7-degrees-21-22.md)。按首项饱和列全 540 个早端形状，使用正系数迫位、整数范数和精确有理区间，最后重算 768 个二元字。
- [任意次数尾部唯一性](../../../proofs/062-c7-saturated-tail-rigidity.md)。给定最大次数、较长侧以及两因子到该次数减 8 的前缀，可行因子对至多一个，其全部系数域也由该前缀生成。

Python 3.13.18。把 verify_degrees_21_22.py 和前两份依赖脚本复制到同一个已创建的 data 盘临时目录，使用启用断言的 Python 标准库执行。

    mkdir -p /Volumes/data/CodexTemp/T2/a30-degree22-check
    cp verify_degrees_21_22.py ../060-gap7/verify_degree_19.py ../060-gap7/verify_degrees_17_18.py /Volumes/data/CodexTemp/T2/a30-degree22-check/
    python3 -B /Volumes/data/CodexTemp/T2/a30-degree22-check/verify_degrees_21_22.py --out /Volumes/data/CodexTemp/T2/a30-degree22-check/new.json

任意次数尾部证书使用 SymPy 1.14.0，从六次式及两侧低端前缀推导权重，并枚举模 163 与 1279 的二元字。

    python3 -B verify_saturated_tail.py --out /Volumes/data/CodexTemp/T2/a30-saturated-tail-check.json

同目录 JSON 为主进程实际复算的输出。
