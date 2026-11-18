# 第66轮八次实二次母式的完整证书

本包排除非有理支持 \(\{2,3,4,5\}\) 及反转 \(\{3,4,5,6\}\)。[中文定理与有限覆盖证明](theorem-zh.md)。

两个实现从88个完整整环域、11,080个 \((L,a)\) 核心重新计算，逐项匹配421个母式和根筛的首拒阶段。其中101个根筛余项的正、反向整数余因子树都到达空层，最晚为空在第92层、第28层。

## 重放

只需启用断言的 Python 3.13 标准库。选择包外尚不存在的输出路径。

```sh
/opt/homebrew/bin/python3.13 /path/to/package/replay.py \
  --out /Volumes/data/CodexTemp/T2/round066-s2345-replay --workers 8
```

重放入口核对清单散列并从零生成六份账本，逐字节匹配冻结结果。输出目录保存各阶段 stdout、stderr、耗时、解释器版本、退出码和输出散列。没有随机参数或浮点容差。

| 文件 | 含义 |
|---|---|
| `source/s2345_entry.py` | 迹根式完整入口、双端相容卷积和精确根筛 |
| `source/verify_s2345_entry.py` | 完整整基独立入口、反向直接反解和精确首拒 |
| `source/s2345_trees.py` | 正、反方向逐层二次域与整数 \(H(GQ)\) 树 |
| `source/verify_s2345_trees.py` | 完整整基的独立双向树 |
| `source/verify_s2345_formulas.py` | 可选 SymPy 1.14 符号复核，独立于标准库重放 |
| `frozen/entry.rows.jsonl` | 421母式及逐项根筛阶段 |
| `frozen/independent-entry.json` | 逐母式首正轴失败点或 Schur 首失败次数 |
| `frozen/trees.rows.jsonl` | 两方向每层活态、拒绝数和前缀散列 |

母式记录为 `[D,L,Tr(a),r,beta,Tr(c),s,Tr(d),t,e,f]`，其中 \(b=La+\beta\)，\(r,s,t\) 为 \(a,c,d\) 在完整整基 \(\{1,\omega\}\) 中的第二坐标。正轴首拒的 `value=[m,n]` 表示清分母后的 \(m+n\omega\)，`side=0` 为 \(P\)，`side=1` 为 \(\bar P\)。整数范数系数按升幂排列。

证书覆盖任意有限次数及任意符号的整数余因子。一般均匀和猜想仍开放。
