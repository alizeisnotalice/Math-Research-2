# 第15轮实际计算与证据记录

日期：2026-11-17（Asia/Shanghai）。本轮一般定理由[证明015](../../proofs/015-selector-old-direction-addback-zh.md)承担；计算用于验证辅助结构、核对严格例及检查真实构造，不替代证明。

## 本轮实际运行

| 执行端与任务 | 实际范围 | 结果 | 时间与并行度 |
| --- | --- | --- | --- |
| 本机非含单位核穷举 | n≤3全部16个合法B/E对，4361个线性子空间F，其中559个乘法闭F；每个配最大F模 | 全部σ(wB+zE)=κ(B) | 0.011715秒，单进程 |
| 本机三角图模型 | 3个z-active双点旧原子、1个inactive双点原子，共16384个取向/标签配置；256个E不稳定C | 指定F分割下界均未超过2κ(C)，最高差−7 | 12进程，计算1.577843秒、驱动2.061897秒 |
| MacBook对旧结果作新比较 | 读取R14的29085个已完成T结果；本轮新算1106个不同B的精确容量 | 26565个等号、1680个差1、840个差2，零违反 | 15进程，1.874993秒 |
| 本机严格实例复验 | gap1/gap2各一个，C、B各遍历全部877个集合分割；T完整DP及两份分割证书 | (κC,κB,κT)=(6,4,9)、(8,6,12) | 0.020127秒，串行，全部exit 0 |
| MacBook原生Metal与CPU | Q=4旧方向加回的3个真实六变量码，各32768个相位 | GPU全Walsh幅值8；98304个相位全部64系数CPU int64一致 | 全部命令缓冲completed，前台进程成功结束 |

本机Python3.15.0；MacBook Python3.9.6、Swift6.4.2、macOS27.3、Apple M5 Pro。所有数学判定均为GF(2)秩、整数动态规划或整数FWHT，零浮点容差，无随机抽样。R15 GPU对照是全量；R14的抽样方式不是本轮范围。

三角图的计算仅对一个指定新代数F算容量下界：每单元局部维数用限制到补集的秩核公式计算，不是T的完整容量上界。一般上界现在已由证明015解决，不继续扩大这一阴性搜索。

## 资源和退出的实际边界

本机12个逻辑CPU配置12个工作进程。监测驱动用独立进程组只统计该批自身进程，一秒一次；采到最高324.6% CPU、RSS求和峰值386784 KiB。RSS求和可能重复计算共享页，不是独占物理内存。任务只有约两秒，不能据此宣称12核持续满载。`triangle-resource.json` 保留实际命令、PID组采样、结束时间及子进程exit 0。没有人为延长负载，也没有停止其他项目。

MacBook的15进程核心容量计算不到两秒，没有本轮持续满载采样。GPU报告的69.623661ms是主机提交、等待、回读区间之和，不是GPU硬件时间戳、内核总时间或占用率。8相位校准为GPU2.863167ms、CPU0.009125ms，不作加速宣传。输出分块缓冲524288字节；没有测全进程GPU峰值内存。

CPU批次完成后工作进程退出；Metal是前台程序，包装器在其成功后写出完整哈希表，完成查询没有残留本批GPU/编译进程。因SSH短暂断开，回执传输重试，但没有重启已完成任务。详见 [GPU回执](RUN_RECORD_Q4_METAL.md)和[远端CPU记录](RUN_RECORD_CORE_CAPACITY.md)。

## 三个Q=4真实取等例

旧变量x=bit0、y=bit2，C=span{1,x,y}，E=B=span{1,a}，a分别为x、y、x+y。新变量为z₀=bit1、w₀=bit3、z₁=bit4、w₁=bit5；基相位为xy+z₀w₀+z₁w₁。

生成器显式写出4份独立h_w、z₀E、z₁a_{w₀}，再加C的对角方向。三个S₄(B,E)都为14维，加回旧方向后均为15维，并含全部全局仿射方向。主Agent另外用位掩码乘法独立重建每份基，和生成器的逐点函数实现一致；重新生成JSON与GPU输入逐字节相同。

κ(C)=κ(B)=4，由本轮上界κ(T₄)≤4+3·4=16；又dimT₄+1=16，故容量恰为16。这里容量取等来自解析夹逼，没有额外声称运行了三个大容量DP。GPU还实际核验每码整个乘法稳定子恰为{0,1}。

新Metal源码来自第14轮已归档4cd0版本，只调整三份新输入的schema、标签元数据及CPU全量对照。源哈希为
`26fba2c735c09609ad144c6233bfb4cfe4f8d31fcb6fe7b8e1d4d3d61a04e8db`，
输入为
`cec2b5f13927ad3365b9918b0a2c8acefa8b21f87e7aa473b5443a115adb5e11`，
报告为
`c0314d3c4bb64dd44e9eed586d48ef5a8665a5f15eacbd45a2a484507645ee0d`。
`q4-metal-manifest-sha256.txt` 还绑定实际二进制；该二进制本身不入Git。环境变量MPS fallback记录不改变执行事实：本轮是原生Metal，无CPU fallback代码。

## 来源与实际核验的区分

七点T容量是**第14轮历史运行**，本轮没有重跑29085个T DP。本轮读取任务/结果，逐一核对输入和历史退出状态，再对1106个不同B新作精确小域DP。两份历史源SHA为d88aaa30…、5221efeb…，完整哈希见 `remote-core-capacity-summary.json`，并与R14独立审计绑定。新比较流水约3.2MiB仍保留MacBook专用目录 `r015/core-capacity-audit/out/per-task-comparison.jsonl`，哈希12768379…；只归档摘要、生成与审计脚本及必要严格例。

两个严格实例另由本机本轮完整复验。C/B的上界算法是定义域全部集合分割枚举，独立于远端码字状态DP；T复用已归档R13精确DP源码，新运行D8、D10两个实例；新旧最优分割均用GF(2)局部秩检查。它们是抽象函数码实例，不被宣称为真实bent反例。

首个远端后处理尝试因Python3.9不支持int.bit_count而失败；替换为bin计数后重跑成功，最终源和报告一致。初期一个子Agent的512例三角脚本遗漏新变量复制，局部支撑筛选条件也相反，所谓正差被主Agent源码审查否定。该输出没有成为候选或数学证据，未归档失败稿；主Agent独立重写后才得到这里的16384例。最终数学论证不依赖这项扫描。

为保留实际执行源码的字节哈希，Swift文件末尾的一行空白原样归档；暂存空白检查仅对本次命令关闭blank-at-eof诊断，其他空白错误仍检查，不修改全局Git配置。

`verify_receipts.py` 最后检查源码/输入/报告绑定、全部计数和Q4独立基重建；它不重新执行GPU或历史大批T DP。两位Agent分别审查非含单位核证明和一般Q递归，最终又补清“继续迭代的删减必须仍含单位”这一量词。

## 复现命令

从本目录运行，使用已有Python；临时数据放data专用目录。远端使用已有 `$HOME/CodexRemoteMath/cache`，不回退到系统盘临时目录。

```sh
PY=/Users/zzh/.local/bin/python3
OUT=/Volumes/data/CodexTemp/T2/a28/replay015
mkdir -p "$OUT" /Volumes/data/CodexTemp/cache/clang
export TMPDIR="$OUT" TMP="$OUT" TEMP="$OUT" PYTHONDONTWRITEBYTECODE=1
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang
"$PY" verify_nonunital_kernel.py --out "$OUT/kernel.json"
"$PY" triangle_graph_scan.py --workers 12 --out "$OUT/triangle.json"
clang++ -O3 -std=c++17 ../round013/capacity_dp.cpp -o "$OUT/capacity_dp"
"$PY" verify_strict_examples.py --inputs strict-examples.json --binary "$OUT/capacity_dp" --out "$OUT/strict.json"
"$PY" generate_q4_input.py "$OUT/q4-input.json"
"$PY" verify_receipts.py
```

MacBook后处理的完整命令见 `RUN_RECORD_CORE_CAPACITY.md`；它依赖R14原始任务/结果，若文件不再存在，可先按R14生成与运行步骤复现，不能仅拿摘要冒充全量输入。最小严格例和Q4输入已独立归档，不依赖原始大流水。

GPU复现：将本轮 `MetalWalshVerifier_r015.swift`、`q4-input.json`、`run_q4_metal.sh` 复制到MacBook的 `$HOME/CodexRemoteMath/a28-selector/r015/q4-metal` 后运行包装器；先保留已存在的历史回执，避免覆盖。Swift编译使用共享模块缓存。运行完成再核对manifest的源、输入、二进制、报告，而不把仅启动命令算作完成。

本轮没有验证目标备份盘内容；提交推送也不代表该备份验收成功。
