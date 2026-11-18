# 第 64 轮本机 MPS 数值记录

本机 MPS 饱和参数化完成 \((32,39)\) 和 \((40,48)\) 两批。每批 12,288 条并行轨迹、36,000 步，随机种子 2026111950、2026111951。使用[第 62 轮参数化源](../062-gap7-numeric/search_gap7_mps_saturated.py)。主进程逐项精确重算全部 576 个保存样本，记录的最小最大残差分别约 0.20482、0.21516。两批均保存并已退出。

    python3 -B ../060-gap7-numeric/verify_saved_samples.py --out /Volumes/data/CodexTemp/T2/a30-mps64-new.json m32n39 m40n48
