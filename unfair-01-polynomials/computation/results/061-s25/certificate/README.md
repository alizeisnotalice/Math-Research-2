# 第61轮八次非有理支持 \(\{2,5\}\) 的排除证书

[中文定理及完整覆盖](theorem-zh.md)。

本包只需 Python 标准库。两个实现从完整实二次整环域重新生成正反前缀，检查全部必要 Newton 约束及严格根环。共81个域、680个 \(a\) 核心、46个通过相容前缀的状态和32个不同母式。根筛只留10项，其正向与反向整数树都实际到达空层。

```sh
/opt/homebrew/bin/python3.13 /path/to/package/replay.py \
  --out /Volumes/data/CodexTemp/T2/round061-s25-replay --workers 4
```

选择包外尚不存在的输出目录，启用 Python 断言。运行会检验清单的散列、从头生成账本并与冻结字节核对。全部阶段的 stdout、stderr、退出码、用时与解释器信息写入新目录。第二套实现扫描完整整基的平移及 \(e\) 整数盒，再以不同方向的 Schur 递降和逐层整数树独立核查。

| 文件 | 含义 |
|---|---|
| `source/s25_entry.py` | 使用高端关系定位整数平移，构建完整母式及精确根筛 |
| `source/verify_s25_entry.py` | 完整整基扫描，独立首拒见证 |
| `source/s25_trees.py` | 正反两方向逐层二次域／HGQ 整数树 |
| `source/verify_s25_trees.py` | 完整整基的独立双树重放 |
| `frozen/entry.rows.jsonl` | 逐母式根筛阶段 |
| `frozen/independent-entry.json` | 首个正轴失败点、嵌入与整基值，或 Schur 失败次数与范数 |
| `frozen/trees.rows.jsonl` | 树活态、拒绝计数与每层散列 |

母式记录为 `[D,Tr(a),r,c,d,Tr(b),s,e,f]`，\(r,s\) 是两个非有理系数在整基 \(\{1,\omega\}\) 的第二坐标。独立正轴见证的 `value=[u,v]` 表示清分母后的 \(u+v\omega\)，`side=0` 表示 \(P\)，`side=1` 表示 \(\bar P\)。范数多项式用升幂整数系数。支撑 \(\{3,6\}\) 由反转同时排除，一般均匀和猜想仍开放。
