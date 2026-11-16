# 第20轮：GPU 终端精确逆纤维

主SMC核以第20轮根矩源为基线；新增独立Metal library，仅在一峰 tau=N+1 时计算规范最短路逆纤维m2、m3。
源SHA-256为 `1f593425af5fe53f98cc3b01109d6498519ad2d555c2ba14332d715c07f26b48`。不改变原随机数、count/CDF/prop/MH流程。

`verification.json` 是n3 P17 L0/4/16和n8 P257 L4的逐group整数/旧数学字段比对。主Agent另从归档源编译、实际运行n8 P257，再以独立根矩trace checker核查前缀/RNG，并与Python逆纤维DP逐组比较，见 `root-independent-replay.json`。n16 P1 onepeak另检查了N=256索引边界。全3×3真实逆像逐词认证见相邻 `r020-inverse-fibre` 及 `r020-cpu-inverse-fibre`。

`trace-compaction.json` 记录五个真实trace仅压缩JSON空白，解析内容完全相同；最大文件约4.45MB。`host_inverse_from_trace.py` 读取这些trace并独立计算m2/m3组和。全前缀合法性由相邻 `r020-gpu-root-moments/verify_root_trace.py` 校验。

n7 P4096 R32 L16只归档一次完整批，`n7-root-recomputed.json` 是主AgentDecimal80重算。R2≈19.48212、R3≈101.31092，对应精确值19.32638、101.64752；跨rep权重ESS≈31.91/32。delta-SE及10%历史诊断门槛均不是置信认证，亦非混合证明。GPU Z是onepeak概率，乘N!才是A1计数；CPU版本直接估计计数。

复现时所有新输出放临时目录（在本目录启动）：

```sh
TASK_TMP=/Volumes/data/CodexTemp/T2/a022-r020-gpu-inverse-replay
mkdir -p "$TASK_TMP"
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules inverse_fibre_gpu_r020.swift -o "$TASK_TMP/inverse_gpu"
"$TASK_TMP/inverse_gpu" 8 257 1 65 65 "$TASK_TMP/n8.jsonl" 1f2e3d4c strat 4 trace
python3 ../r020-gpu-root-moments/verify_root_trace.py "$TASK_TMP/n8.jsonl.trace.json"
python3 host_inverse_from_trace.py "$TASK_TMP/n8.jsonl.trace.json"
python3 summarize_inverse_gpu.py n7-p4096-r32-t50-l16.jsonl --output "$TASK_TMP/n7-summary.json"
```

n12正式批另在临时目录运行；未完成数据不在本轮结果之内。
