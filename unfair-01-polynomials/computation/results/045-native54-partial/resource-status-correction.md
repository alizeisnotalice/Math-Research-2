# 资源状态更正（2026-11-17）

此前“无自有本地/远端计算进程残留”的表述超出了可验证证据，特此撤回。

- 已证实：本地 degree-54 queue controller `32751` 已退出；本地 task18
  scheduler/controller 也已退出。
- 未证实：雷电桥物理断连且严格 SSH 不可达时，无法查询 task18 的远端
  native 进程、CPU controller `75942` 或其 workers。因此这些远端进程是否仍
  存活、是否占用 GPU/CPU 或内存，均为 **unknown pending reconnection**。
- task18 只具有本地 `partial_error` receipt（0 verified batches）；这不是远端
  已终止、也不是内存已释放的证据。task19--23 未由本地 queue 登记为完成。

恢复严格连接后，必须先只读核验这些远端 PID、journal/final 和资源状态；在此
之前不得重跑或将 task18 记为 terminal/complete。
