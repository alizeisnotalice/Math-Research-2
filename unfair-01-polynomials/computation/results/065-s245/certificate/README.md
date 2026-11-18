# 第65轮八次非有理支持 \(\{2,4,5\}\) 的排除证书

[中文定理及完整覆盖](theorem-zh.md)。

本包只需 Python 标准库。两个实现从81个完整整环域、680个 \(a\) 核心重算通常与例外两类尾部，得到262个完整母式和逐项相同的根筛阶段。35个根筛余项的正、反向整数树均出现空层。

```sh
/opt/homebrew/bin/python3.13 /path/to/package/replay.py \
  --out /Volumes/data/CodexTemp/T2/round065-s245-replay --workers 4
```

选择包外尚不存在的输出路径，并启用断言。脚本验证清单散列，从零生成各账本，再与冻结结果逐字节对照。各阶段的 stdout、stderr、退出码、耗时和解释器信息写入输出目录。

| 文件 | 含义 |
|---|---|
| `source/s245_entry.py` | 迹根式的完整入口、两端相容卷积与精确根筛 |
| `source/verify_s245_entry.py` | 完整整基的独立尾部、Newton 与首拒见证 |
| `source/s245_trees.py` | 正反方向逐层二次域／整数 \(H(GQ)\) 树 |
| `source/verify_s245_trees.py` | 完整整基的独立双向树复算 |
| `frozen/entry.rows.jsonl` | 每个完整母式的根筛阶段 |
| `frozen/independent-entry.json` | 逐母式首正轴失败点或 Schur 首失败次数 |
| `frozen/trees.rows.jsonl` | 正反向每层活态、拒绝数及散列 |

母式记录为 `[D,Tr(a),r,b,Tr(c),s,Tr(d),t,e,f]`，\(r,s,t\) 是三个非有理系数在完整整基 \(\{1,\omega\}\) 中的第二坐标。独立正轴见证的 `value=[m,n]` 表示清分母后的 \(m+n\omega\)，`side=0` 表示 \(P\)，`side=1` 表示 \(\bar P\)。整数范数系数按升幂排列。

定理同时排除反转支持 \(\{3,4,6\}\)，覆盖任意有限次数的整数余因子。一般均匀和猜想仍开放。
