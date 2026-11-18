# 第59轮八次非有理支持 \(\{1,3\}\) 的排除证书

[中文定理和完整覆盖](theorem-zh.md)。

本包无需数据下载或第三方模块。两个 Python 标准库实现分别用半整数根式坐标和完整整环坐标。它们从所有允许的实二次域重新生成前缀，直到第4行的共同整数母系数约束。独立实现产生866条准许两侧单位区间的共同平移，每条给 \(C_4=2\)，所以完整候选数为0。额外根筛和尾部树没有输入。

```sh
/opt/homebrew/bin/python3.13 /path/to/package/replay.py \
  --out /Volumes/data/CodexTemp/T2/round059-s13-replay --workers 4
```

选择包外不存在的输出目录。Python 须启用断言。复跑将从头求两套账本，与冻结结果逐字节对照，并保存两个进程的退出码、标准输出、标准错误和耗时。

| 文件 | 含义 |
|---|---|
| `source/s13_exact_four_rows.py` | 以根式坐标生成前缀，先施加 \(C_4=1\)，再检验共同 \(d\) |
| `source/verify_s13_integral_basis.py` | 以完整整环基重算，先构造共同单位平移，再直接算 \(C_4\) |
| `frozen/four-row-entry.json` | 主运行各必要入口计数 |
| `frozen/four-row-entry.rows.jsonl` | 完整母式的空账本 |
| `frozen/independent-entry.rows.jsonl` | 866个共同 \(d\) 的精确见证 |
| `frozen/independent-entry.json` | 独立计数和见证哈希 |
| `manifest.json` | 源及冻结账本 SHA-256 |

见证每行是 `[[D,Tr(a),r,c,Tr(b),s,d],C4,A4,B4]`。系数 \(a,b,A_4,B_4\) 使用 \(\mathcal O_K=\mathbb Z+\omega\mathbb Z\) 的整基，\(A_4,B_4\) 以 \([u,v]\) 表示 \(u+v\omega\)，\(r,s\) 分别是 \(a,b\) 的 \(\omega\) 坐标。记录中的共同 \(d\) 未使用第四阶 Newton 界，也未限制到系数盒。

本定理处理任意有限整数余因子；非有理位置 \(\{5,7\}\) 由反转一并排除。一般均匀和猜想仍开放。
