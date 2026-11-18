# 次数 23 至 26 的缺口排除

[次数 23、24](../../../proofs/063-c7-degrees-23-24.md) 覆盖 8,640 个早端分支。只在次序正确的仿射子系统使用 Cramer 法，五个残余二次式用模 163 互素证书排除。

[次数 25、26](../../../proofs/063-c7-degrees-25-26.md) 覆盖 138,240 个早端分支，最后在两个模数核验 3,061 个晚端字。冻结的 396 个显式有限域多项式恒等式可只读核验，在过渡到剩余域扩张时仍有效。

在本目录执行以下命令，用启用断言的 Python 3.13.18。前一部分为标准库，后一部分使用 SymPy 1.14.0。

    mkdir -p /Volumes/data/CodexTemp/T2/a30-degree26-check
    cp ../060-gap7/verify_degree_19.py ../060-gap7/verify_degrees_17_18.py ../062-gap7/verify_degrees_21_22.py \
       verify_degrees_23_24.py verify_degrees_25_26.py degrees-25-26.json /Volumes/data/CodexTemp/T2/a30-degree26-check/
    cd /Volumes/data/CodexTemp/T2/a30-degree26-check
    python3 -B verify_degrees_23_24.py --out degrees-23-24-new.json
    python3 -B verify_degrees_25_26.py --witness degrees-25-26.json --out degrees-25-26-new.json

两份主进程重放记录分别是 degrees-23-24.json 和 parent-identity-replay.json。程序可省略 --witness 重新生成 25、26 次的乘子。
