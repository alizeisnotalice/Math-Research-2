# r012：原条件律按第二出生秩 \(\tau\) 的整数 SMC

这是辅助数值核，不是远峰概率问题的证明，也不提供构造或反构造。
`DESIGN.md` 固定了原图、均匀下降标号、\(E_\tau=\{T_2=\tau,K=2\}\)
和未归一化估计对象 \(Z_\tau,U_\tau\)。其条件期望归纳只用**每个父粒子**的
条件期望 offspring 数（等价地，所有 child 平均的条件均值）；不要求每个 child
同分布或粒子独立，也不把伪随机数写成数学随机性。
最后一秩不重采样。跨 \(\tau\) 的 \(\sum U/\sum Z\) 只是数值比值。

## 实际核验

- `original_tau_smc_audit.py` 穷举 n=2,3 的均匀排列与受限树。得到
  \(A_2(2)=8,Q_2=1/2\)，\(A_2(3)=143808,Q_3=1100/2247\)，所有非零
  \(\tau\) 层逐项一致。该文件的完整输出是 `tau_smc_audit_final.json`。
- n=2 GPU smoke 精确输出 \(Z=1/3,U=1/6\)。n=3，P=512、每层64重复的
  合并比值为 0.4934823093；这是有限重复的数值诊断，非精确值。
- n=7，P=16384、每层8重复、\(\tau=2,\ldots,49\)：
  \(\sum\bar Z=2.2621706428138447\times10^{-10}\)，精确
  \(A_2(7)/49!=2.2535264802956476\times10^{-10}\)，相对差 +0.00383584；
  合并比值 0.3451438143，精确 \(Q_7=0.344418314818\)，差 +0.0007254995。
  这是固定种子有限 PRNG 模拟的一次锚定，不能充当误差界；同种子重跑只核验
  可复现性，并不构成新的独立重复。
- n=8，P=32768、每层16重复、\(\tau=2,\ldots,64\)：
  \(\sum\bar Z=9.603664010052252\times10^{-15}\)，合并比值
  0.3369305513。此包没有 n=8 的精确 \(A_2,Q\) 基准，故不声称准确度。
  \(\tau=63,64\) 各16次全灭，均作为零估计保留；非零层中最小
  **单步不同直接父索引数** 是26（\(\tau=62\)）。它不是跨代祖先数。

完整 n=8 运行的 `run_end.wall_seconds` 为 212.987011 秒；所有重复的
Metal GPU 时间戳和为 66.175797 秒，同步 command-buffer 提交加等待和为
136.744382 秒。后二者是跨重复累加，不能与全进程墙钟简单相加或相除；
GPU 时间戳/墙钟的覆盖比也不是整机 AGX 利用率，不能据此声称硬件被持续占满。
n=7 在未插桩版本运行，未报告硬件 GPU 时间，不能倒推。

## 计重与日志更正

每一步令 \(a_i=|A_r(X_i)|\)，每个 child 从 \(0,\ldots,\sum_i a_i-1\)
均匀抽整数。每256粒子 GPU 精确求组和；CPU 仅从组和定位组；GPU 在该组扫描。
区间长度给出精确的边缘父选择概率 \(a_i/\sum_j a_j\)，顶点选择使用64位拒绝法。
child 的 mutation 流由 `(replicate, rank, child-index, domain)` 的新键产生，
不会复制父的未来 RNG 状态；这并非对现实 PRNG 独立性的证明。

早期 n=7 emitter 将死粒子输出为 `logZ:-inf`，这不是 JSON。原始文件
`data/n7-p16384-r8.jsonl` 原样保存；`normalize_dead_logz.py` 仅在记录同时
满足 `dead=true,Z=0,U=0,min...=0` 时把该单个 token 改为 `null`，否则拒绝。
它转换16行，所得合法文件与固定种子重跑逐字节相同。所有旧 raw 的
`min_distinct_ancestors` 字段实际含义只是最小单步不同直接父索引数；校正摘要
改用 `min_one_step_distinct_parents`。

## 完整标号调试及边界

生产 n=8 State 只有已访问掩码和两根，不能恢复完整标号，故不从其输出声称
独立峰检验。`history_debug_audit.py` 是**单独 CPU 模拟**，并非实际 GPU 全标签
轨迹审计；它保存 n=4、P=64 的完整下降序列；
它保存26条可完成轨迹（\(\tau=2,\ldots,14\) 每层两条），每条由独立的终态
严格局部峰检查器验证 \(K=2\)、出生秩为 \(1,\tau\)、并核对根距。\(\tau=15,16\)
在每层128次小池尝试中无完成轨迹，记录为有限搜索失败，不是不可行性证明。

`source/tau_smc_integer_guarded_future.swift` 把可表示的掩码范围显式限制为
\(2\le n\le12\)，并检查 P、R 和 \(\tau\) 参数；该未来版本未替换已冻结生产。

## 复现

在 Apple Metal/Swift 环境（本轮设备为 Apple M6）中：

```sh
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules
swiftc -O -framework Metal source/tau_smc_integer_n8_timed_used.swift -o tau_smc
./tau_smc 8 32768 16 2 64 n8.jsonl
python3 source/summarize_tau.py n8.jsonl > n8-summary.jsonl
python3 source/normalize_dead_logz.py data/n7-p16384-r8.jsonl n7-legal.jsonl
python3 original_tau_smc_audit.py
python3 source/history_debug_audit.py history-debug.json
```

运行后应按 `HASHES.sha256` 核对冻结源、原始数据和摘要。所有计算程序已退出；本包不含二进制和缓存。
