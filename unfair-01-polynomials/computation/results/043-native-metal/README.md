# 第43轮 MacBook 原生Metal接管

旧v3严格归属的控制器和Wolfram内核已退出；exit为intentional partial/rc255，不伪装完整。v1/v2/v3已收到的连续并集715,915,264 masks，经父独立区间和字节哈希核验。

新native从task13的11,272,192起点接续，至33,554,432，共22,282,240 masks/170批，首批与余批均有complete回执。父从原始roots+C直接重构保存的全部680个候选，均有数值负系数，最大乘积目标误差1.17404e-6；容差1e-5。只重构保存候选，不声称证明每个mask的结果，不作严格排除。

源码Swift/Metal为97918305c6c2bcc02db02c36a84bf970264a5d387c2491cf3cae4f7233a3711d；实测remote binary为8235c31db552c331e3552044347d3ba784029e7ab7f027bb2440cff23b0ea457。首批controller总耗时19.077秒、余下169批24.767秒，含SSH/传输/编译等开销，不与纯GPU kernel基准混淆。Apple M5 Pro环境及编译器版本沿用42轮已实测工具记录。

每任务独立input/source/binary与回执绑定；builder重算真实已收区间，无缺口或重叠才能接续。缺失coverage flag、错误hash/count等负测已拒绝；早期首批旧格式保持原样，派生range另绑定原ledger和receipt。未完成任务仍按partial/unknown处理。

复现候选独验：在本目录执行 `python check_native_task13.py`。原始inputs.json.gz由该脚本直接读取。
公开包的三个range-extract只保存实际区间事件投影，其原始完整流SHA单独登记，不冒充原流字节副本。完整旧流留在本地原路径；对完整旧流运行audit_handoff.py可复现父哈希/区间核验。公开投影可重算连续区间及计数，但不能单靠投影复核原流字节。
正式调度器保留实际执行源码，其INPUTS和SSH路径需按已授权设备配置；不要运行旧覆盖范围来冒充新结果。

后继task14控制器在SSH超时后自行退出为partial_unverified_stream，远端journal尚未取回，范围未知；14..23未纳入本轮完成范围。待查远端真实状态后接续，不重复未知范围。CPU51/52和本机49独验也单独运行。
