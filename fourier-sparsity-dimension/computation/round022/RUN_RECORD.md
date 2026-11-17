# R22 实际运行与复现

本机Python 3.15.0；Apple M6、macOS27.2、Swift6.4.2。MacBook Apple M5 Pro、macOS27.3、Swift6.4.2。全部整数或Fraction精确计算，零容差，无随机种子。

- type2_core.py：完整16384种14点核心符号，恰128合法，组成以22为代表的唯一仿射陪集。
- gauss_and_cross_audit.py：mu3..6的二次型补集Walsh、谱支撑上限、Gram严格正定下界，以及第二型参考信号与整数能量等号条件。
- cross_partition_audit.py：mu3..6共5440个完整内在平移，逐项检查dyadic谱片覆盖，并用F2消元计算参数法向秩。此为局部模型，非全局bent分割实例。
- mm_cross_witness.py：六组(mu,r)=(3,0),(3,1),(3,2),(4,0),(4,1),(5,0)，直接核验两个真实相位的全部整数Walsh及差支撑重量。
- 两机相同Swift源码各检查16384掩码、1048576有符号系数。CPU使用定义式14项点积，GPU使用六阶段int32 FWHT；全部零分歧。128合法掩码完整集合与Python一致。

源文件哈希、实际命令、设备、计时及编译/运行退出码见各JSON与两份receipt。主Agent实际运行verify_receipts.py交叉验收。远端首次程序因将仿射陪集误比为线性空间而exit10，最终源码修正该验收条件后完整重跑exit0；初次运行不作为最终证据。782336个目标集合外的系数是预期候选否决。

本机编译显式指定/Volumes/data/CodexTemp/cache/swift-module-cache。远端实际编译命令未显式传模块缓存，不能声称已传；回执保留实际命令。每批4096相位，两个主要GPU缓冲各1MiB。结果取回后远端自有进程已退出；本机程序同步执行结束，未打开MATLAB，无持续满载声明。

复验存档先运行`python3 verify_receipts.py`。各Python源单独运行会更新同目录JSON的计时及内容；应在副本复现以保留归档回执。Swift复现示例：

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-module-cache R22Type2CoreMetal.swift -o /Volumes/data/CodexTemp/T2/a28/r022-local/reproduce
/Volumes/data/CodexTemp/T2/a28/r022-local/reproduce /Volumes/data/CodexTemp/T2/a28/r022-local/reproduce.json
```

全维分类、刚性与预算由解析证明承担，有限检查只验证相应实现及实例。
