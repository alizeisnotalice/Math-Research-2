# 固定第 31 项的本机 MPS 探索

PyTorch 2.15.0 MPS 对 \((39,33)\) 完成 12,288 条轨迹、36,000 步，固定第 12 至 31 项的优化目标。保存的 288 组 Float32 系数以二进制有理数重算卷积。运行已保存并退出。

    python -B ../060-gap7-numeric/verify_saved_samples.py --out /Volumes/data/CodexTemp/T2/a30-mps39n33-check.json m39n33

本批为数值探索。较低次数的严格结论来自独立的有限域覆盖。
