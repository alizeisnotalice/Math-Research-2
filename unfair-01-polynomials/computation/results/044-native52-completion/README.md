# degree52数值任务最终覆盖与父独验

24个任务的已收区间总和1,308,622,848 masks，全部连续到各自终点；这只是数值排名覆盖，不是严格排除。前v1/v2/v3和task13已在前轮记录，本轮封存task14..23。

父独立重算原始日志并集、每任务实际manifest/input/source/binary/receipt、逐批真实journal范围及保存mask绑定，再从原roots+C重构全部保存样本。源Swift固定97918305…、binary8235c31…，并非只接受恢复器成功标记。实际耗时10.926秒。

本轮父Python3.13.18在product tolerance1e-5下结果：{"nearunknown": 113, "negative": 17295}。unknown全部保留；没有近似数转成精确反例或排除。
代理旧Python3.15回执在部分mask得到不同残差和分类；同一输入/同一checker/hash/mask列表用3.13重跑已复现父结果。旧3.15记录只作不一致诊断，不称完全独立一致。父数值记录是本轮最终展示依据；浮点版本差异更不支持严格结论。

正常完成与SSH断线后的恢复分别保存原始与recovered receipt，缺失旧receipt不伪造补写。实际文件hash、覆盖区间和远端启动都另绑定；旧恢复器缺少的coveragefalse/预期binary等防伪门槛由父在本次独验明确检查，后继degree54另版本加固，不改写已执行源。

公开包保存task14..23所需原始输入、manifest、journal/final与ledger/receipt；大件逐文件gzip并有archive-bindings。原始24目标inputs和独立roots重构器沿用../043-native-metal；完整前v1/v2/v3原流保留在本地，公开范围投影不冒充其原字节。parent checker保留真实本机路径，复现需按索引恢复相应输入或调整路径参数。

本批所有本题native内核已确认退出；新degree54为不同目标、独立source/inputs/队列，未计入本轮覆盖。
