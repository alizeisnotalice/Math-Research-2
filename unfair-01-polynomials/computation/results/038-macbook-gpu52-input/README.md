# MacBook52次后继GPU输入与实际启动

24个52次目标，seed2026111616；目标0为精确公平控制。GPUArray Real32、L64、batch32768。实际以每目标独立Wolfram15内核 -run 参数传入约217KB程序，remote_gpu.py逐个启动并wait回收，整批14400秒到时保存已返回前缀并回收当前内核。远端仅RAM，结果stdout回本机data；本轮仍运行，不宣称完整排除。

两次长stdin入口未进入搜索；TERM未能结束首个内核后定向KILL，PID26240与28213已确认不在进程表、两个bridge均rc255已回收。保留两个失败bridge和exit元数据，失败不计完成。单目标-run的32768mask小试实际成功后才启动完整批次。公开脚本不包含SSH密钥/配置。

实际五次ioreg GPU Device Utilization采样74,0,73,100,99；记录为gpu-telemetry.json。瞬时达到100%不等于全程无空档；下一优化将减少CPU评分传输，并探索GPU驻留迭代。运行未结束的样本留独立temp目录，不以快照替代完整输出。
