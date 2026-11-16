# 全GPU有序前缀移动核与n7校准

soa_frozen.swift为实验n≤16实现：264-byte compact状态，加两个device ushort[P*256]前缀缓存。count/CDF/选父/复制真实父前缀/添点/MH/终端矩均在GPU上；host每秩同步用于错误检查和可选trace。死亡后move早退，终端无重采样和移动。

主Agent从冻结源独立编译并实际重跑62例（n3/n8/n16，P17及257/513跨组），root-frozen-62-audit.json保存准确输入与原输出哈希。verify_soa_frozen.py以独立Python集合、整数RNG和原始birth定义重建每个活粒子，核验真实GPU全部32+32状态字、父索引、分组CDF、每步移动及终端整数矩。logZ/Z/U/Q按非零纯相对容差5e-12核对；死亡0与一峰null分开。

归档只保留7份代表trace，其JSON空白被压缩；原/归档哈希映射及解析对象完全相等检查见trace-normalization.json。其余原始大trace位于临时卷，完整复现程序包含62例。

```sh
python3 reproduce_all_62.py
python3 verify_soa_frozen.py n16-p17-t2-l4.jsonl.trace.json
python3 summarize_gpu_move.py n7-p4096-r32-l0.jsonl n7-p4096-r32-l16.jsonl
```

完整n7 P4096,R32,τ2..49，L0/L16同域8200000000000000；每批1536记录均完成并退出。L16的Q≈.344408667，精确Q≈.344418315；分层delta SE≈.000476735。L0/L16同seed是配对诊断，不能当独立批次合并。核审计和有限校准不认证混合速度或渐近概率。

排除的旧审计错误已修正：host伪补高字、terminal未核、onepeak内部dummy矩误当零、n32入口超256前缀容量、向量zip截断、极小Z/U使用绝对容差。最终核验重新冻结，并实际拒绝截断parent/word以及篡改Z/U；旧证书未冒充最终版本。

复现时请先确认 data 挂载，将 TMPDIR/TMP/TEMP 及编译输出设到 `/Volumes/data/CodexTemp/T2`；Swift/Clang 模块缓存为 `/Volumes/data/CodexTemp/cache`。环境：Python 3.15，Apple Clang 21，Swift 6.4.2/Metal。源码相对路径调整不改变数学核；原运行源定位保留在本轮报告，归档哈希见 manifest。有限核验不外推到所有 n。

本轮另完成 n16 P16384、L16、R128 的 τ2 与 τ257 两个单层试验，domain 分别为8400000000000000与8900000000000000。独立重复权重 ESS 分别13.2738、5.6280，最大单次份额23.12%、27.70%；仍有强权重集中，不能将两批之比当作已认证的 A₂/A₁。τ2 的加权 Q≈.284356 仅属于该出生秩，不能代替完整条件律。两批分别约321.8/318.5秒，均正常退出、零swap。`n16-pilot-summary.jsonl`由 `../r018-million-pilot/summary_logz_deadzero.py` 对完整128条逐文件计算；该脚本核对seed、死亡和覆盖。旧百万粒子完整原核批次继续运行，不纳入本轮完成结果。
