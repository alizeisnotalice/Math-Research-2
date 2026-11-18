# 第 63 轮双机系数探索

本机 MPS 的次数对 \((28,31)\) 运行 16,384 条并行轨迹、36,000 步，种子 2026111813。保存样本 288 行逐项有理重算，最小记录残差约 0.21398，进程已释放后退出。源沿用[第 62 轮饱和参数化](../062-gap7-numeric/search_gap7_mps_saturated.py)。

MacBook Apple M5 Pro 的原生 Metal 版本实现同一低端形状和两首项所迫零位。使用 Swift 6.4.2，先真实执行 64 条轨迹、256 步的试跑，再执行 \((35,39)\)、16,384 条、32,000 步，种子 2026111816。运行约 19.32 秒后正常退出。两次各保存 16 行，主进程在以 IEEE Float 重新量化后逐项计算精确有理卷积。更宽的 35 对次数扫描仍在运行。

    python3 -B ../060-gap7-numeric/verify_saved_samples.py --out /Volumes/data/CodexTemp/T2/a30-m28n31-new.json m28n31
    python3 -B verify_metal_samples.py smoke.json --out /Volumes/data/CodexTemp/T2/a30-metal-smoke-new.json
    python3 -B verify_metal_samples.py production35-39.json --out /Volumes/data/CodexTemp/T2/a30-metal-production-new.json

保存的舍入样本均未达到原方程的精确等式。
