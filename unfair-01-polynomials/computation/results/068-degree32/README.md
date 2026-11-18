# 第 32 行内部卷积的完整精确证书

[定理与各片覆盖](../../../proofs/068-c7-degree32-full.md)。

Python 3.13.18、SymPy 1.14.0。重放从压缩档恢复有向区间账本和 16 片单位见证，重新生成原卷积式，核查 27,196 条有限域单位恒等式与 1 条实数正下界。主进程的实际回执为 parent-verification.json。

    python -B replay.py --out /Volumes/data/CodexTemp/T2/a30-j32-replay --workers 6

临时输出须放在包外新目录。原文件归档散列、区间及程序绑定见 degree32-archive/verification.json。两次调度重切均追加新片后重放，原来完整的前段保留。
