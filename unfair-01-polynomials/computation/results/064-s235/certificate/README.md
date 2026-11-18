# 第64轮八次非有理支持 \(\{2,3,5\}\) 的排除证书

[中文定理及完整覆盖](theorem-zh.md)。

本包只需 Python 标准库。从88个完整整环域重新生成5,440个 \((L,a)\) 核心，遍历通常与例外尾部分支。两套实现分别得到相同的91个完整母式及必要根筛阶段。16个根筛余项的正、反向整数树均实际到达空层。

```sh
/opt/homebrew/bin/python3.13 /path/to/package/replay.py \
  --out /Volumes/data/CodexTemp/T2/round064-s235-replay --workers 4
```

选择包外尚不存在的输出路径并启用断言。脚本检查包内清单散列，从零生成全部账本并与冻结结果逐字节对照。所有阶段的 stdout、stderr、退出码、用时及解释器信息写入输出目录。

| 文件 | 含义 |
|---|---|
| `source/s235_entry.py` | 迹根式上下端精确覆盖与根筛 |
| `source/verify_s235_entry.py` | 完整整基独立构造、尾部扫描、首拒见证 |
| `source/s235_trees.py` | 正反方向的二次域／整数 \(H(GQ)\) 复核 |
| `source/verify_s235_trees.py` | 完整整基的独立逐层复算 |
| `frozen/entry.rows.jsonl` | 完整母式及必要根筛阶段 |
| `frozen/independent-entry.json` | 独立计数及第一失败点评估或 Schur 首失败次数 |
| `frozen/trees.rows.jsonl` | 两向活态、拒绝数与整前缀散列 |

每个母式记录为 `[D,L,Tr(a),r,u,c,Tr(d),t,e,f]`。\(r,t\) 是 \(a,d\) 在完整整基 \(1,\omega\) 的第二坐标，\(b=La+u\)。独立正轴见证中 `value=[m,n]` 表示清分母后的 \(m+n\omega\)，`side=0` 表示 \(P\)，`side=1` 表示 \(\bar P\)，整数范数系数按升幂排列。

定理覆盖反转支持 \(\{3,5,6\}\) 和任意有限整数余因子。一般均匀和猜想仍开放。
