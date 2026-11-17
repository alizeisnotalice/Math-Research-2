# MacBook 52 次根组 GPU 续跑 v2

实际 Mathematica GPUArray 小批试跑通过，内核退出已检查；生产从已审计旧批的边界续接，严格 SSH 使用原主机密钥。新收到1372批、44,957,696个掩码后 SSH 返回255，没有终止完成事件，远端内核是否退出尚待重连核实。此批只记部分数值探索；保存行还需独立重建。状态见 ../041-compute-snapshot.json。

原始已生产脚本的 B 逆变换公式正确；错误仅在从未启动的旧续跑准备稿。v2 修正后，主线程实际运行 audit_inverse_formula.py，对六个已知掩码用 CPU 多项式乘积独立重建，最大误差约1.85e-5，小于3e-4门槛。具体参数、hash和输出见 parent-formula-regression.json。没有把浮点误差门槛作为精确反例证书。

input.json.gz 为续跑输入；make_resume_plan.py 构建边界，remote_gpu_resume.py 与 search_resume_v2.wl 执行任务。bridge_v2.py 的 all_owned_processes_reaped 仅指本地 SSH 进程，遇255不证明远端进程已经退出。当前不得盲目再次启动同一输入。

可在本包目录用 Python 运行 audit_inverse_formula.py 检查公式；完整调用参数见脚本及随附回执。已退出的试跑与正在核实的生产分别记账。Swift/Metal 优化准备尚未作为本包生产结果。
