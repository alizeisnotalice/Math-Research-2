2026-11-18，Python 3.15.0，整数 Walsh 与 GF(2) 枚举。归档目录复现命令：

    python3 min_support_4d.py ../round035/exhaust_4d_local.py min_support_4d_result.json
    python3 verify_e4_extension.py e4-extension-audit.json
    python3 enumerate_mm_e4.py mm_e4_enumeration_result.json

四变量标准基相位的维数1至4取等码数为60、120、30、1500。枚举先考察合法二次投影，再枚举投影子空间、仿射核以及全部线性抬升，并独立核验取等见证。另对固定四原子组件穷尽65536种额外方向，其中128个相容，8个在组件内，其余120个外重均至少4。

本机及MacBook独立枚举一个指定六变量MM模型的5760个允许置换与相位，得到1474560个组合，模组件后有90个拓展达到最小外重8。MacBook Apple M5 Pro 的CPU与Metal对这些拓展各16个参数的全部92160个系数逐项比对，零分歧。编译、运行均exit0，含编译用时约1.85秒，进程已退出。模型内计数与解析全维下界分别记录。复现Swift时输入 MacBook 枚举结果与输出路径，使用已确认的模块缓存。完整终端回执保存在 data 临时目录 r036-mac。

第二重量间隙沿用第9轮 proofs/009b-residual-budget-zh.md 的定理与来源。一般穿孔距离和四维下界的依据见本轮证明。
