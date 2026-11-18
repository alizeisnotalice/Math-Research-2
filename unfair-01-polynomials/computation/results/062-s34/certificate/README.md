# 第62轮八次两位置 \(\{3,4\}\) 的排除证书

[中文定理与完整覆盖](theorem-zh.md)。

本包只需 Python 标准库。两个实现从完整实二次整环重新生成正反前缀，计算必要 Newton 界和严格根环。扫描287个平方自由底数、保留55,526个 \((L,M,a)\) 核心、35,546个不同母式。根筛留下2,454项，其前向树均实际出现空层，最大首次空层53。反向树有7项仍活到深度上限，均由前向证书排除。

```sh
/opt/homebrew/bin/python3.13 /path/to/package/replay.py \
  --out /Volumes/data/CodexTemp/T2/round062-s34-replay --workers 4
```

选择包外尚不存在的输出目录，启用 Python 断言。运行会检验清单的散列、从头生成全部账本并与冻结字节比较。各阶段保存实际 stdout、stderr、退出码、用时和解释器信息。两套实现的共用数学约束见定理。独立实现使用不同整环坐标、完整整基取整及相反方向的 Schur 递降。

| 文件 | 含义 |
|---|---|
| `source/s34_entry.py` | 主精确入口、两端相容卷积及根筛 |
| `source/verify_s34_entry.py` | 完整整基的独立前端、尾部及根筛复算 |
| `source/s34_trees.py` | 正向和反向卷积／HGQ 树 |
| `source/verify_s34_trees.py` | 整基独立重建所有树层 |
| `frozen/entry.rows.jsonl` | 逐母式首拒阶段 |
| `frozen/trees.rows.jsonl` | 逐母式树的活态、拒绝数与每层散列 |
| `frozen/independent-entry.json` | 完整母式一致性回执 |
| `frozen/independent-trees.json` | 完整树一致性回执 |

逐母式记录 `[D,L,M,Tr(a),r,u,c,d,e]`，\(r\) 是 \(a\) 在完整整基 \(\{1,\omega\}\) 的第二坐标，第四次母系数 \(b=La+u\)。`root_pass` 只是必要根筛保留。树回执中的方向未决与原模型状态分离，正向空层已覆盖每个母式。

反转后的支持 \(\{4,5\}\) 同时排除，一般均匀和猜想仍开放。
