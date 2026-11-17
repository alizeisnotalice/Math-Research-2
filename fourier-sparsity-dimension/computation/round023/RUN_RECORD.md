# 第23轮实际运行

本机Python3.15.0、macOS27.2、Apple M6、Swift6.4.2；MacBook macOS27.3、Apple M5 Pro、Swift6.4.2。整数/F₂运算，零容差。

## 两片真实MM族

参数(μ,r)=(3,0),(3,1),(4,0)，分别256、1024、1024相位。Python实际检查每个相位全Walsh为±q、部分Walsh为0或±q及选片包含关系，并精确消元求法向投影及核限制秩。独立审查员也实际重跑了最终Python源。

两机相同最终Swift源分别用CPU定义式字符和与GPU整数FWHT逐项比较全谱及部分谱，共每机1081344有符号系数，零分歧。两机全部结果行完全一致，与Python结构字段相同。两份receipt记录实际命令、源/输出/二进制哈希及退出码，均exit0；两机编译均显式传共享模块缓存路径。本机执行结束，远端取回后核实自有进程退出，未打开MATLAB。

曾修正Swift拼写及把q/2误当谱幅度的验收条件：每片点数为q/2，但非零谱幅度为q。最终源在两机完整重跑，旧失败输出未归档。各计时见JSON与回执。

## 完整spread构造

主Agent在归档路径实际重跑最终desarguesian_cross_partition_audit.py。k3检查全部2048相位；k4检查512个块码相位及5120个仿射基平移，共5632次（含重复），不是全部17维码；k5检查4097个固定种子20261117的码内样本，不是全部27维码。生成秩依次11、17、27；有限域模多项式、块重量、代表块及谱值均记录于JSON。最终源明确验证所采样前d个生成元独立。全维成立性由符号和恒等式证明，不依赖有限抽样。

## 复现

先执行`python3 verify_receipts.py`核验存档绑定。Python程序会覆写同目录JSON计时，应在副本复现以保留归档。Swift示例：

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-module-cache R23TwoPatchMetal.swift -o /Volumes/data/CodexTemp/T2/a28/r023-local/reproduce
/Volumes/data/CodexTemp/T2/a28/r023-local/reproduce /Volumes/data/CodexTemp/T2/a28/r023-local/reproduce.json
```

有限核验用于检查实现与实例；不把GPU一致性或哈希一致性写成一般证明。不存在本轮持续满载声明。
