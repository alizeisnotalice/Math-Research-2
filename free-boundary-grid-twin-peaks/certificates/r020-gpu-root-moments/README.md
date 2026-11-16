# r020：终端根向量矩扩展审计

基线是只读文件
`../r019-gpu-move/soa_frozen.swift`（实际归档路径见第19轮清单）
（SHA-256 `a8b234543d556f7152020772a81c63f8a50341dbd94af7abb53e0ce87387fc56`）。
本目录的 `root_moments_r020.swift` 只改 terminal kernel 和 terminal-output
路径：每一 256-particle group 额外写出

\[
F_0=\sum_i a_i(\|t_{p_i}\|^2+\|t_{q_i}\|^2),\qquad
G_0=\sum_i a_i\langle t_{p_i},t_{q_i}\rangle,
\]

其中 (t_v=2(x_v,y_v)-(n-1,n-1))。`F_0` 用 `UInt64`，`G_0` 用 `Int64`。
它不读取随机数；count/CDF/prop/move 的 Metal source、buffer indexes 和 dispatch
次序未改。`r020-vs-frozen.diff` 是完整差异记录。

对双峰（`tau<=N`）每个 replicate 输出

\[
D=Z D_0/T_N,\quad F=Z F_0/T_N,\quad G=Z G_0/T_N,
\]

以及 `F_over_Z`、`G_over_Z`、`V=F/(8Z)`、`C=G/(4Z)`、
`Q=D/(NZ)`（既有 `U_over_Z`）。独立整数核验严格检查
`F_0-2G_0=4D_0`；浮点输出的对应残差只允许机器舍入量级。`tau=N+1` 时所有这些
根矩字段为 null，且 GPU root group buffers 被置零。

## 实际审计

`r020-audit-certificate.json`（5/5 成功）记录精确命令和每个 trace SHA-256：

- n=3, tau=2, P=17，L=0、4、16：冻结 r019 与 r020 的 `seed_hex,Z,U,logZ,U_over_Z,scaled_U,dead`
  JSON numeric token 完全相同；r020 全 trace 的父代、prefix、state、terminal D/F/G
  由独立 Python 重算通过；
- n=16, tau=2, P=257, L=4：跨两个 workgroup 的 GPU CDF 与 root F/G group trace 独立通过；
- n=3, tau=10=N+1, P=17, L=4：一峰输出中的 U/Q/D/F/G/归一化根矩字段均为 null，
  trace 独立通过。

随后 `r020-frozen-62-certificate.json` 以同一 source/verifier/binary SHA-256
实际覆盖 r019 的全部 57 个核心例及 5 个跨 workgroup 例；62/62 的 GPU run 和
独立 verifier 都成功，且每一记录包含完整命令、output SHA-256 与 trace SHA-256。
`terminal-integer-bounds.json` 给出 n=16 下单 group 的解析上界
`F0<=230400`、`|G0|<=115200`，故 UInt64/Int64 均无溢出风险。
`root-numeric-negative-tests.json` 记录 F/G 置零、翻号、加倍六项修改全部被拒绝。

## n=7 回归校准

`r020-n7-calibration.json` 是以 `n=7,P=4096,R=32,tau=2..49,L=16,domain=8200000000000000`
实际运行的 1536 条 replicate 汇总。它与 r019 归档逐 token 比较
`seed_hex,Z,U,logZ,U_over_Z,scaled_U,dead`，得到 0 个不匹配。
按先对每个 tau 的 32 条取均值、再取比的结果为

- `V/n² = 0.11260576092626597…`（相对 r018 精确值误差 `-4.1761e-5`）；
- `C/n² = -0.05959857254742299…`（误差 `-3.6937e-5`）；
- `Q = 0.34440866694737793…`（误差 `-9.6479e-6`）。

相应分层 delta-SE 仅作描述性诊断：V/n² `1.0593e-4`、C/n²
`1.6356e-4`、Q `4.76735e-4`；它们不是认证置信区间，也不证明混合。
所有非死亡 replicate 的浮点 `F-2G-4D` 相对残差最大约 `2.44e-16`。

入口 guard `P²N<=UInt64.max` 下，最坏 n=16 时 `P<=268435455`；所有 P 的
terminal host 累加仍满足 `F0<=241591909500`、`|G0|<=120795954750`，故也安全于
UInt64/Int64，而不只是单 group 安全。

复现（在本目录启动；输出和编译放临时盘）：

```sh
TASK_TMP=/Volumes/data/CodexTemp/T2/a022-r020-root-replay
mkdir -p "$TASK_TMP"
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules \
  swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
  root_moments_r020.swift -o "$TASK_TMP/root_moments"
python3 replay_trace_catalog.py --binary "$TASK_TMP/root_moments" --output-dir "$TASK_TMP/traces"
python3 summarize_gpu_roots.py --help
python3 check_terminal_integer_bounds.py
```

`replay_trace_catalog.py` 从既有62例证书读取实际命令参数，重新运行GPU并独立核对；旧时间记录不作为新运行。
`summarize_gpu_roots.py` 验证完整分层覆盖、种子、死亡记录与矩恒等式，再汇总比例和分层delta误差；能从logZ恢复浮点下溢的非死亡权重。
`summarizer-checks.json` 记录截断、重复、错种子、错矩和错死亡标记的拒绝测试，以及合成下溢输入的比例不变性检查。

当前 SHA-256：

- `root_moments_r020.swift`: `0ddb392f0c9b28a8a40e78944b71f77e3ee0e591a95d63b3a0b2ffbe692addeb`
- `verify_root_trace.py`: `851deadfb3955f6a41a37208d817c242f9bdeb8f24c41715028b4d1dc91e6211`
- `r020-audit-certificate.json`: `c8b182f9d15cac7df79c2acd150394adf68791f8fada7ea773c2da3bb3180688`
- `r020-frozen-62-certificate.json`: `04e4875f84b73bfb226939fe67ec09a7af814ea59fb3a01e2998b5f21034ed20`
- `r020-n7-calibration.json`: `ead6fedf7af705cfcaed0744b33f0425e5b24f009ce82423a01811de3afe2acc`

该审计验证实现与终端矩记账，不是关于 MH 混合速度或双峰条件概率的认证。
