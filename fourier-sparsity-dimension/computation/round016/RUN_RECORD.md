# 第16轮实际计算与证据

日期：2026-11-17（Asia/Shanghai）。本轮结论见[证明016](../../proofs/016-full-selector-extension-rigidity-zh.md)。以下区分：小维真值表的完整枚举、由结构定理得到的高维扩张分类、以及对预测构造族的实际CPU/GPU全量核验。

## 本机种子与独立容量检查

`build_extension_seed.py` 用与字符真值表的精确Hamming计数计算Walsh系数，完整检查全部65536个四变量相位，得到896个bent相位。固定基准相位x₀y₀+x₁y₁后，允许的相对方向对全部32个仿射方向封闭；模五维RM(1,4)恰有28类。

程序从零子空间开始，逐个加入允许商方向并检查整个新陪集，完整枚举允许的商子空间。各维数计数1、27、45，三维层为空；所有包含极大者恰为45个七维扩张。每份128码字均重新检查bent性，支撑分割DP给容量8。

`verify_base_capacities.py` 使用第二种算法：每个七维码中完整枚举2825个含单位线性子空间，筛乘法闭代数并求其最大稳定模块。45份合计127125个子空间、4335个代数，容量全部8，与DP一致。本机使用12个进程，最后归档运行耗时0.261131秒。没有持续CPU利用率采样，不称12核满载。

种子程序还对二变量最大仿射码的四变量满选择器，完整检查全部65536个任意外加方向，恰128个共同合法，且正好等于原七维选择器。边界例中，先删成常数分支再恢复全局仿射时，xz可合法加入，但不与满选择器的yw共同合法。另一个七维例RM(1,4)+span{xz,yz}的全部128相位、256个三变量分支均被检查：每支恰四个等幅频点、二次部分秩2、两支频谱不交覆盖。8条参数频面公式也逐项正确。

最后一次种子报告运行耗时0.221434秒，版本Python3.15.0。开发中为加入严格归约及容量交叉检查而重跑过本机小程序；归档的是最后源码对应的报告，不把这些重复计为更多独立输入。`extension-seed.json` 的896相位内容在此期间保持字节不变，与两项远端任务严格绑定。

## 六变量预测族的实际CPU/GPU运行

最小输入 `extension-seed.json` 包含排序的896个16位相对相位掩码，SHA-256：
`c94010a703187ddeccc15f94d3390e7777324cd2d99742a1476872a0b61f3df0`。
对每个有序三元组(a,h₀,h₁)，其中a∈{0,1}，构造

```
g = h0 | (h0 << 16) | (h1 << 32) | (h1 << 48)
if a == 1: g ^= 0xffff0000ffff0000
F0 mask = 0x935f6ca06ca06ca0
```

这与位顺序x₀,x₁,y₀,y₁,z,w上的相位h_w(x)+az一致。该族共有1,605,632个不同方向，**整个族不是线性空间**。独立分支的全仿射平移闭包使该族对十一维选择器的所有生成元封闭。

| 执行 | 本轮实际范围 | 结果与时间 |
| --- | --- | --- |
| MacBook C++ CPU | 15线程、全部1605632相位、每相位64点整数FWHT | processed=expected；坏相位/坏系数/最大幅值偏差全部0；内核0.166073秒、驱动1.194532秒；exit 0 |
| MacBook Swift Metal | 同一JSON、784块×2048相位、每相位64点int32 FWHT | 全部命令缓冲completed，坏系数0；主机提交/等待/回读区间之和1367.080665ms |
| Metal启动校准 | 首2048相位的全部64系数，与本进程CPU int64逐项比较 | 全等；GPU主机区间4.207584ms、CPU1.207375ms |

两个完整程序独立验证同一族的全部Walsh幅值。**没有声称全族CPU/GPU有符号系数逐项比较**；该更强比较仅对首2048相位实际执行。Metal滚动值`8c8eb544649d0383`由源中offset1469598103934665603定义，不是标准FNV初值，也不是加密摘要或独立全系数证书。

CPU最小文本输入、JSON到文本的提取回执、C++源码、运行器、结果及软件版本齐备；`cpu-extension-receipt.json` 包含实际返回码。其结果哈希为262be721…，详见[CPU记录](RUN_RECORD_EXTENSION_CPU.md)。GPU源码与报告经 `metal-manifest-sha256.txt` 绑定，环境、完整分块记录及前台退出见[Metal记录](RUN_RECORD_EXTENSION_METAL.md)。

MacBook环境为Python3.9.6、Apple clang21.0.0、Swift6.4.2、macOS27.3、Apple M5 Pro。GPU执行是原生Metal，未使用PyTorch/MPS；环境中的fallback变量不是后端证明。所有数学判定均使用精确整数或GF(2)，随机种子不适用，零浮点判定容差。FWHT中间绝对值≤64，不存在int32溢出。

## 资源、完成与传输

远端CPU任务短于0.5秒采样周期，只保留开始/结束样本，没有持续满载证据。Metal的计时不是GPU硬件时间戳、内核时间或利用率，不能用来宣称GPU满载或加速。每块输出缓冲为524288字节，未测整个进程峰值内存。没有为维持负载重复已完成任务。

CPU首个SSH握手超时发生在发出运行命令前；查询确认没有已启动任务后才真正运行一次。两个远端目录互相独立，CPU和GPU都成功结束，完成查询无本轮计算/编译进程。本机小程序与12进程审计也均exit 0。没有停止其他课题，百万相位没有保存成大型本机原始数组。

## 结果依赖与范围

- 896个旧相位、45个旧包含极大扩张及其容量8，是本轮小维完整计算的结论。
- 指定十一维六变量子码的1,605,632个共同合法方向、2025个包含极大扩张、全部十五维且容量16，由**结构定理加小维完整分类**得到。没有穷举2^64个六变量真值表，也没有将此结论用于不含指定子码的其他六变量码。
- CPU/GPU全族扫描实际认证预测族的bent性；任意未列方向的排除依赖解析证明。45个最大扩张容量还由与DP不同的全部代数枚举复验。
- 两位Agent独立审查满分支刚性、整个扩张空间的投影论证、异质容量公式及只保留常数时的四谱归约。参数—频率联合支撑非仿射的明确例子由本机重新计算。一般A28仍开放。

`verify_receipts.py` 检查最终源码/输入/输出绑定、种子数据、完整CPU计数、全部784个GPU分块及报告范围；它本身不重跑百万相位。备份盘副本没有核验。

## 复现

从本目录执行，使用已有Python和编译器。大输出及编译缓存必须位于data临时目录；MacBook使用既有专用远程目录和共享缓存。

```sh
PY=/Users/zzh/.local/bin/python3
OUT=/Volumes/data/CodexTemp/T2/a28/replay016
mkdir -p "$OUT" /Volumes/data/CodexTemp/cache/clang
export TMPDIR="$OUT" TMP="$OUT" TEMP="$OUT" PYTHONDONTWRITEBYTECODE=1
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang
"$PY" build_extension_seed.py --out "$OUT/extension-seed.json" --report "$OUT/seed-verification.json"
"$PY" verify_base_capacities.py --input "$OUT/seed-verification.json" --out "$OUT/base-capacity.json" --workers 12
"$PY" extract_seed_masks.py --input "$OUT/extension-seed.json" --out "$OUT/seed_masks.txt" --receipt "$OUT/seed-extraction.json"
clang++ -O3 -std=c++17 -pthread verify_extension_cpu.cpp -o "$OUT/verify_extension_cpu"
"$OUT/verify_extension_cpu" "$OUT/seed_masks.txt" 12 "$OUT/cpu-result.json"
"$PY" verify_receipts.py
```

MacBook CPU实际使用15线程。若复现完整资源回执，把相同文件复制到一个新的远程计算目录、编译，再运行 `run_extension_cpu.py`；它使用自身目录和15线程、等待子进程并记录返回码。先保留历史回执，勿重放覆盖。

GPU复现：将 `ExtensionMetal.swift`、`extension-seed.json`、`run_extension_metal.sh` 放入MacBook的 `$HOME/CodexRemoteMath/a28-selector/r016/extension-metal`，运行包装器。它使用共享Swift模块缓存；完成后核对manifest及全部分块记录。不要仅凭启动命令声称完成。
