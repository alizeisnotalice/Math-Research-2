# 第 31 行内部卷积全部支撑排除

[定理、逐行约束与复算命令](../../../proofs/067-c7-degree31-full.md)。

把归档中每片的 gzip 部件按编号解压拼接。从 source 目录用 Python 3.13.18、SymPy 1.14.0 重放，旧依赖已经在相邻的 round065-gap7-next 目录按原始散列冻结。输出应使用包外的临时目录，移植时请保持两个源码目录同层。

主进程重放 27,648 个分支的有向区间，按四片执行 13,214 条有限域单位见证，并重新推导唯一双模余项的实平方恒等式。结果见 parent-verification.json。完整压缩账本及各文件 SHA-256 见 degree31-archive/verification.json。

包内一次重放入口使用全新输出路径，自动核对各压缩部件及全部源码散列。

    python -B replay.py --out /Volumes/data/CodexTemp/T2/a30-j31-replay --workers 4
