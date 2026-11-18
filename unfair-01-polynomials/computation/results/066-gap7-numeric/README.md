# 第 7 项缺口的本机 MPS 样本

本机 PyTorch 2.15.0 MPS 完成 \((33,31)\) 的 12,288 条轨迹和 36,000 步。保存的 288 组 Float32 系数用二进制有理数重算。该批保存后已退出。64 次排除由另行完成的精确证明承担。

    python -B ../060-gap7-numeric/verify_saved_samples.py --out /Volumes/data/CodexTemp/T2/a30-mps31-new.json .
