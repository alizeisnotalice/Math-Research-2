# 第21轮实际运行记录

本机 Python 3.15.0、macOS 27.2 arm64；原生 Metal Apple M6。MacBook 为 Apple M5 Pro、macOS 27.3、Swift 6.4.2。全部使用整数运算，容差为零，无随机输入。

主 Agent 归档验收时实际重跑四个 Python 程序，并由归档源重新编译运行本机 Metal；本机编译指定共享模块缓存。逐次耗时见结果 JSON，本机编译和执行回执见 local-receipt.json。远端实际最终执行及源码/输出绑定见 macbook-receipt.json；该次远端编译未显式指定模块缓存。两机程序都结束，无持续满载声明；未打开 MATLAB。

| 检查 | 实际覆盖与结果 |
| --- | --- |
| 六点符号 | 64种，32种合法，恰一个仿射限制陪集 |
| 四分之一幅度网格 | 9^6=531441，恰零向量与32种饱和符号 |
| 一个自由变量 | 全4096种符号，64种局部平台信号，均有共同自由字符 |
| 两机 Metal | 每机531441向量、8503056有符号系数，与直接16×6 CPU点积零分歧 |
| 真实6+6+4分割 | D=7全部128个bent相位；七柱体逐点等于部分谱，法向/参数投影秩均4 |

Metal 采用四阶段 int32 FWHT，每批4096向量；两个主要GPU缓冲合计524288字节。系数绝对值至多24，无溢出。7939296个不在目标集合的系数属于预期候选否决，不是GPU错误。初始Swift多行字符串闭合及通过条件已修正，归档源实际执行通过。两机33个解的完整集合与Python独立枚举一致。

复现：在本目录运行 Python 程序（各自输出同目录JSON），再编译 Swift 源并传入输出路径。重新运行会更新计时；若复验存档回执，先运行 `python3 verify_receipts.py`，避免覆盖其输出哈希。例如本机编译：

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-module-cache R21LocalMetal.swift -o /Volumes/data/CodexTemp/T2/a28/r021-metal-local/reproduce
/Volumes/data/CodexTemp/T2/a28/r021-metal-local/reproduce /Volumes/data/CodexTemp/T2/a28/r021-metal-local/reproduce.json
```

六点幅度刚性与任意自由维数的结论由解析卷积/Parseval证明承担；这些有限运行核对实现和具体实例，不替代一般证明。
