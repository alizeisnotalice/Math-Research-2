# 第14轮实际计算记录

日期：2026-11-17，Asia/Shanghai。一般 A28 与一般旧方向加回容量不等式均未解决。
本目录保留最终源码、最小输入、摘要和必要证书；不把数值阴性结果当作一般证明。

## 执行范围与实际结果

1. 本机先完成四点预检：全部16个含单位函数码、64个合法代数对，精确容量均不超过两倍基码容量，1.44519秒。每份最优分割另以完整码字与 GF(2) 秩检查。
2. MacBook 枚举七点集上全部2825个含单位函数码、52897个合法代数对。23812个代数稳定整个码，适用证明014的包含关系；其余29085对取最大稳定子码 B，完整运行精确 DP，全部成功、零超时、零遗漏、零反例。15个并行子进程，批次190.818秒；最大提升维数16。对任意更小 B 的覆盖来自容量单调性。重复映射到64个坐标是满射，保留函数代数及容量，不是概率权重归一化。
3. 固定四变量 D=7 三角 bent 码的全部115个含单位代数：每个取最大模块，作 Q=2 选择器再加回全局仿射方向。最终六变量码维数分布为9:8、10:24、12:6、13:24、14:39、15:14。MacBook 15路精确 DP 耗时1.17188秒，容量分布10:8、11:24、16:83。这里并非115个码全部维数15；不依赖不规范基表示的全局去重断言。
4. 同一115份输入在 MacBook **原生 Swift Metal GPU** 上完整检查1,347,584个相位，所有整数 Walsh 系数绝对值恰为8。CPU int64 对照33,680个相位的全部64个系数，逐项一致；另外8相位校准也一致。每码整个乘法稳定子均仅有0与1两个元素。GPU完整遍历是本次实际执行，不是只检查设备可用。
5. 用户要求本机同步参与后，新增三个八点模型：基码是前四点偶重、后四点自由的七维码，κ=12。其三种代数、最大模块及提升基由 `eight-point-candidates.json` 完整给定。所有码都得到κ=24。本机 D19 使用12线程，内核报告59.3514秒；MacBook 两个 D18 各15线程，内核报告16.9326/15.3183秒，驱动含监测的墙钟为19.777/17.270秒。远端先复验历史 D7、D13 模型，容量8/16，均 exit 0。
6. 三个八点结果均有独立 GF(2) 分割证书和直接和上界核验：本机8+4+4+4+4；远端 free-pair 为8+8+4+4；mixed 为12+4+4+4。程序逐项检查分块投影稳定整个码，再对各小分量分别运行旧串行与新并行 DP，两者一致。这给出三个指定实例的独立上界，不是一般加回定理。
7. MM 单仿射方向容量跳跃例在 q=4、8 时分别为6→8、10→16。先在 MacBook Python3.9.6实际执行；最终归档检查发现该脚本把换行写成字面 `\\n`，导致输出文件不是合法 JSON。修复序列化后，本机 Python3.15.0重新实际执行，0.013559秒，结果一致。归档 `mm-affine-jump.json` 是**修复后的本机运行**，不是先前远端运行。

所有运算使用 GF(2) 秩、整数 DP 或整数 FWHT，零浮点判定容差。搜索为确定性穷举；GPU与CPU对照首码全查，其余每码以固定LCG种子 i+20261117（零起算码序号i，再按源码取 seed|1）选8个相位，非随机研究输入。浮点数仅用于计时与资源采样。

## GPU 来源、计时与边界

设备：Apple M5 Pro；系统 macOS27.3 (26C5035r)，Swift6.4.2。实际源
`MetalWalshVerifier.swift` SHA-256 为
`4cd0f699f659608d369a2d18c4ed2fb06d07f0025c36e34a8c784c170a4ce87b`；
输入为 `selector-input.json`，SHA-256
`c092419514114cf44d9e5e7d85c2d6d4dc998483d02e6f55c325bce1fdf1609c`；
报告 SHA-256 为
`7b66b4145da87b8fd81a5ba7958e7678e9ff39a18fbba3919d7589c741b2e313`。
这三者及二进制的绑定见 `metal-bundle-sha256.txt`、`cpu-gpu-binding-audit.json`。

每个64线程组做一个64点 FWHT；每块2048个相位，输出缓冲524288字节。
中间值绝对值至多64，不存在 int32 溢出。所有命令缓冲实际达到 completed；源码没有 CPU fallback 路径。
`PYTORCH_ENABLE_MPS_FALLBACK=0` 只是环境记录，该执行器不用 PyTorch/MPS。
`total_gpu_ms=993.121547` 是主机 commit、等待和回读区间之和，**不是 GPU 硬件时间戳、占用率或完整进程耗时**。
8相位微型校准为 GPU 2.255875ms、CPU 0.006166ms；不能据此宣称 GPU 加速。
未测GPU占用率或进程峰值内存。后续 Agent 曾改动一处错误元数据并另报计时；该源/报告未形成此处完整匹配包，因此不混入证据。

## CPU、内存与退出

- MacBook 是主要存储和计算端，Python3.9.6、Apple clang21.0.0、15逻辑 CPU。七点批次76次资源采样中，整机最低 idle=0%，本批子进程采样峰值1034.7%。另一次在15个 DP 子进程活动时实测整机90.13% user、9.86% sys、0% idle。这证明该批期间整机达到满载，不能把整机负载全归属于 A28。
- 后续两个 D18 任务都短于60秒采样周期，只采到启动时自身129.3%/146.0%，不声称持续满载。
- 本机 D19 确实启动12线程并完成输出，独立证书通过，PID61091在退出后查询已不存在。原启动器未收集子进程 wait 返回码，故**该次原始 exit code 未观测**；不补写成 exit 0。小分量核验器和其他有 wait 的子进程均实际返回0。
- 本机16:43:24的89.33% user、10.66% sys、0% idle是在该 D19 已结束后的整机采样，不能作为它跑满 CPU 的证据。两机都实际使用，但没有持续同步满载的观测。
- 用户要求释放本机内存后，将大批数据和计算留在 MacBook，只回传小证书。普通 `purge` 被系统拒绝；`sudo -n purge` 需要密码，也未成功。观察到本机可用内存比例21%、29%、62%，不把变化归功于失败的 purge。没有停止其他课题进程。
- MacBook CPU 子进程全部返回0；D18 驱动最后进程查询为空。Swift 前台程序结束，`set -e` 包装器在其成功后才写完整哈希表，没有后台 GPU worker；进程退出释放自身 Metal 资源。本机本轮任务均已结束。

详细记录另见 `resource-and-exit-receipt.json` 和 `RUN_RECORD_D18_REMOTE.md`。

## 完整性、失败修正与归档范围

`verify_batch_results.py` 在 MacBook 上逐条独立检查全部29085条原始记录和最终115条记录：ID唯一完整、输入完全相同、退出码0、JSON结果一致、分割互斥覆盖、摘要哈希一致。它没有独立重算全部局部秩；另用 `verify_selected_partitions.py` 对20条边界见证完整重算局部码维数，全部得κ20。八点见证用低内存的投影秩验证器。

主Agent最终 `verify_receipts.py` 检查归档源码/输入/报告绑定、GPU/CPU相同输入、所有相位数和直和证书摘要。该检查不是重新跑GPU或29,085次DP。

七点完整结果25,166,766字节和任务10,326,300字节留在 MacBook
`/Users/zhengzhihao/CodexRemoteMath/a28-selector/r014/`，分别为 `results-rerun.jsonl`、`tasks.jsonl`；不提交完整原始流水。
其哈希分别为
`5221efeba0c8a4435a88c8bf388d42e31f23a29eb1d5e6e17677fcb5208382e6`、
`d88aaa304d30126a0a1338fa06a10d0a00026abbd68cdca9dd403dba9fd6494c`。
生成器、全覆盖审计、摘要与边界最小输入均归档；115条完整小结果也归档。备份盘副本未核验。

首个七点驱动尝试因相对可执行路径丢失 `./` 在启动 DP 前失败，修复为绝对路径后才有本轮29085个成功记录。
早期选择器构造漏掉两个变量副本，产生非法rank19输入；已丢弃，最终输入由修复生成器在MacBook重建，字节哈希一致，且完整GPU Walsh核验通过。不得将这些失败输入计入覆盖。

## 复现

以下从本目录执行；`PY` 选择已有 Python3.9以上，`OUT` 必须设到本机 data 临时目录或 MacBook 远程专用计算目录。不能在仓库中生成大批原始流水。
本机编译缓存使用 `/Volumes/data/CodexTemp/cache`；MacBook 使用已有 `$HOME/CodexRemoteMath/cache`。

```sh
PY=/usr/bin/python3
OUT=/Volumes/data/CodexTemp/T2/a28/replay014
mkdir -p "$OUT" /Volumes/data/CodexTemp/cache/clang
export TMPDIR="$OUT" TMP="$OUT" TEMP="$OUT"
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang
clang++ -O3 -std=c++17 ../round013/capacity_dp.cpp -o "$OUT/capacity_dp"
clang++ -O3 -std=c++17 -pthread capacity_dp_parallel.cpp -o "$OUT/capacity_dp_parallel"
"$PY" build_addback_tasks.py --points 7 --out "$OUT/tasks.jsonl" --summary "$OUT/generation.json"
"$PY" run_capacity_batch.py --tasks "$OUT/tasks.jsonl" --binary "$OUT/capacity_dp" --results "$OUT/results.jsonl" --summary "$OUT/summary.json" --workers 15 --batch-seconds 600
"$PY" verify_batch_results.py --tasks "$OUT/tasks.jsonl" --results "$OUT/results.jsonl" --summary "$OUT/summary.json" --binary "$OUT/capacity_dp" --out "$OUT/audit.json"
"$PY" verify_selected_partitions.py "$OUT/summary.json" --out "$OUT/boundary.json"
"$PY" build_selector_input.py --out "$OUT/selector-input.json"
"$PY" verify_mm_affine_capacity_jump.py --out "$OUT/mm.json"
"$PY" verify_receipts.py
```

本机/远端并行度分别可选12/15。新编译二进制哈希可能因平台不同而不同，不要求与历史二进制字节一致。
115个任务可直接从归档 `selector-capacity-results.jsonl` 的每个 `task` 字段提取、按 `id` 排序写 JSONL，再交同一批处理器，`--batch-seconds 300`。

GPU复现：将最终Swift源和 `run_metal_walsh.sh` 复制到 MacBook 的
`$HOME/CodexRemoteMath/a28-selector/r014-gpu`，将 `selector-input.json` 复制为该目录的
`selector-input-canonical.json`；运行 `zsh run_metal_walsh.sh`。包装器使用共享Swift缓存，前台完成后写 `manifest-sha256.txt`。
先保留已有完整回执，避免重放覆盖原始证据。

八点复现：对输入中每条记录调用 `capacity_dp_parallel THREADS D <basis整数列表>`，捕获标准输出及返回码。远端 `run_d18_remote.py` 同目录需有源和输入，仅运行两个 D18；其内置900秒上限并等待每个子进程。每个结果再运行：

```sh
"$PY" verify_eight_point_case.py --inputs eight-point-candidates.json --index 0 --result local-eight-point-result.json --serial-binary "$OUT/capacity_dp" --parallel-binary "$OUT/capacity_dp_parallel" --out "$OUT/eight-point-audit.json"
```

另两个依次用 index 1/2 和对应 `remote-eight-point-*-result.json`。整数DP的层间同步确保只读取更低重量的已完成状态；没有跨层未完成读。
