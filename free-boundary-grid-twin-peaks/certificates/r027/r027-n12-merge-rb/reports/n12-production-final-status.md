# n=12 merge-RB r026 生产批：退出与审计状态

状态：**正常完成且输入/实现不变量审计通过；有限诊断，不是混合、SMC 精度或原问题结论认证。**

生产参数为 `n=12,P=131072,R=32,tau=2,L=16,method=strat,domain=9e260c0000000001`。完整 raw 为 metadata + 32 replicate + 32 `firstmerge_rb_calibration` + run_end，共 66 行。驱动、`time` 包装、`caffeinate` 和 GPU 子进程均已消失；telemetry 的末记录为 production PID 57428 以 returncode 0 退出。

先行小 trace (`n=12,P=17,R=1,L=16`) 已由独立 Python 逐 rank 重放 144 ranks，并重新构造末端完整 word、K=2、L/距离/root/group/raw 字段；RB payload 独立审计也通过。它只核对这条有限实现路径，不能证明链已混合或样本估计准确。

生产 raw 使用加固版 `aggregate_production_r026.py`，其 SHA-256 为 `f1593cc93db21d778289966c80539606b079d702fa3fd5c197394c27d4737389`；该审计拒绝重复 ID、缺 metadata/RB、错误 seed、非有限数、组整数恒等式、RB histogram 与 run_end 等错误。输出审计 SHA-256 为 `698412c507b5ba754d97e354e3d899c3643a5cc55d7abc01026e718c6cfc320c`，并报告 `audit_passed: true`。

有限加权诊断：normalizer ESS 30.9294、最大份额 0.0541990；m-lift ESS 30.9500、最大份额 0.0537680；`Q_mRB=0.316566061858082`，旧 L-group 路径 `Q_L_legacy=0.316446842883523`，`A2/B` 的 RB 估计为 2.16358503441652。以上均不得表述为误差认证或一般规律。

资源：wall 119.30 s，max RSS 1,237,450,752 bytes，peak footprint 1,580,500,648 bytes，0 swap、0 received signals。已完成计算的子进程已退出，未遗留本批 GPU 内存占用。

核查材料：raw SHA-256 `05fde7b272cd8eabed5fc62f81371033b7a380178a9a578e4372ae3b913ee7f1`；派生源 SHA-256 `3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6`；二进制 SHA-256 `5de42d7ec2eb946cfb805b59df658cf2f25d87e6018d4561d3e7a7b3a2a9d598`；receipt、telemetry、资源文件哈希均由 receipt 固定。
