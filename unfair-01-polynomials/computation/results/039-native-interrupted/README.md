# 原生GPU梯度批次：仅已收到的进度

m=60,n=100,L=256；8,192并行初值，计划12,000步，种子2026111617。只收到10条进度至1,000步，最后quartic loss为15958.609375；未收到完整系数、最终状态或退出回执。本机重启中断结果流，不能声称远端生产批次完成或据此排除反例。SSH恢复后已确认旧控制器/内核不在。

小规模smoke为已完成独立批次，其回执不替代生产退出。`verify_prefix_gradient.py` 的Float64有限差分误差6.13e-8只核验梯度公式；四项生产进度负测试通过。主线程实际重跑 `audit_native_gradient_partial.py --root . --out <data临时路径>/audit.json`，结果见parent-native-audit.json。不得把进度loss当作精确候选。
