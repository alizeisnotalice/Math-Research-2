# n=7,p1 终止的只读诊断（r028）

诊断范围：仅检查 `PID 91005`、其 `/usr/bin/time` 父进程 `91004`、以及
2026-11-17 04:49--04:53 UTC（12:49--12:53 CST）针对该 PID、可执行名或输出名的
系统日志。未重跑程序、未终止任何进程，也未读取无关应用的日志内容。

## 已确认事实

`receipt.json` 记录 p1 自 04:48:58 UTC 左右运行至 04:52:18 UTC，wall time
200.429 s，以 `exit_code: 1` 结束，且没有产生 `n7-v24-p1.json`。标准错误在完整
写完 layer 1--12（最后一条为 layer 12）后只有：

```
time: command terminated abnormally
time: signal: Invalid argument
```

`/usr/bin/time -l` 对该进程报告：max resident set size 1,949,925,376 B，peak memory
footprint 2,873,657,768 B，**0 swaps**，113.28 user s、141.41 sys s 和 3,558,723
involuntary context switches。作为比较，成功的 p0 完整跑完 layer 44，用 85.03 wall s，
max RSS 1,903,689,728 B，0 swaps。两者均为同一冻结 source/binary：

```
source  5f93fa39695e0f4badabcbabaad2befab920fa1590ea2337fe47a09201587db9
binary  e57308e1d422b949ed43c39f8939f3fb046d7041e6f1adf61a93d48b7994c37c
```

统一日志中唯一明确提及该 PID 的内核记录是 12:49:48 CST：

```
process punctured_layer_[91005] caught waking the CPU 45001 times over ~50 seconds,
averaging 895 wakes / second and violating a limit of 45000 wakes over 300 seconds.
```

这是一条高唤醒率警告；其文本不说该进程被终止，也没有给出信号或内存理由。

## 未发现及结论边界

在上述窄时间窗对 `processID == 91005 OR processID == 91004`，以及对字符串
`91005`、`punctured_layer_mod`、`n7-v24-p1` 的 unified-log 查询中，**未找到**该任务的
jetsam、OOM、crash 或 explicit termination 条目。相应时间窗内的 DiagnosticReports
文件名中也无该任务名。这只能说明已查询来源没有给出原因，不能排除系统层资源压力、
外部信号、运行时/内核异常或未被保留的日志。

因此本轮严格结论为：**p1 的实际终止原因未知；现有记录不构成 OOM 或 jetsam 归因的
证据。** p1 没有 checkpoint，不能作为一个有效模残数。高内核态时间、上下文切换和高
唤醒率是后续重试时应降低并发/减少唤醒的性能线索，而不是已认证原因。

## 原始证据和完整性指针

原始运行材料未归档到本证书包；本目录只记录诊断和下列完整性哈希：

| 文件 | SHA-256 |
| --- | --- |
| receipt.json | 596cccb1c2e939da1a157434a3ac8ca03d409a7cd2e62d8f16e74a8e2dc44355 |
| p0.resource.log | 4b53d4e572558c1f3395088135f669d98f156ebc1f2cb28b9a8e1d2aed3ce100 |
| p1.resource.log | 36d22587fc704f582a25ba5f1dd80fc9333e0a019364a123ce9ca92c97e88a42 |
| n7-v24-p0.json | 241e53cb24fe1d8cd017d49484d8f2ba476f28368897fc40c31ac7d12660cadb |

用于只读复核的命令：

```
/usr/bin/log show --style syslog --start '2026-11-17 04:49:00+0000' \
  --end '2026-11-17 04:53:00+0000' \
  --predicate 'eventMessage CONTAINS[c] "punctured_layer_mod" OR eventMessage CONTAINS[c] "91005"' \
  --info --debug
```
