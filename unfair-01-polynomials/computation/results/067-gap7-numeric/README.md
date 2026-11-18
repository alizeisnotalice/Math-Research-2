# CPU/GPU 混合运行的已归档数值批次

MacBook Apple M5 Pro 原生 Metal 完成 40 对次数。每对 32,768 条轨迹、32,000 步，合计保存 640 组 Float32 样本。本机重建有理卷积并核对全部残差。该批的监督进程已退出。

本机 PyTorch 2.15.0 MPS 完成 \((47,41)\) 的 12,288 条轨迹、36,000 步。优化阶段固定 \(C_{12},\ldots,C_{30}=1\)。对 288 组保留样本作二进制有理数卷积核验，最小保留二元残差约 0.22919。运行已保存并退出。两套样本仍用于数值途径选择。

    python -B verify_metal_batch.py mac-sweep40 --scratch /Volumes/data/CodexTemp/T2/a30-mac40-rerun --out /Volumes/data/CodexTemp/T2/a30-mac40-check.json
    python -B ../060-gap7-numeric/verify_saved_samples.py --out /Volumes/data/CodexTemp/T2/a30-mps47n41-check.json m47n41
