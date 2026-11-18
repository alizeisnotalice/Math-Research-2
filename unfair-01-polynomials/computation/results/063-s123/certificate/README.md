# 第63轮八次非有理支持 \(\{1,2,3\}\) 的排除证书

[中文定理和有限覆盖](theorem-zh.md) · [中心第6行解析引理](analytic-center-zh.md)

中心引理在纯解析层面证明高端 \(s-tg\ne0\) 分支无解。剩下 \(s=tg\) 时，由完整实二次整环和两嵌入 Newton 界构造所有可能的低端前缀。两套整数实现实际考察了1,140个 \(a\) 核心及181,088个带 \(g\) 的两系数前缀。独立实现直接构造出1,199个同时符合两边单位区间约束的公共第四行母系数，每项均给 \(C_4=2\)。合法第四行数为零，因而没有额外根筛输入。

```sh
/opt/homebrew/bin/python3.13 /path/to/package/replay.py \
  --out /Volumes/data/CodexTemp/T2/round063-s123-replay --workers 4
```

需要启用断言，输出目录选择包外尚不存在的路径。标准库代码会检查包内散列，从零重算两个实现并与冻结字节对照。各阶段 stdout、stderr、退出码、用时及解释器信息均记录到输出目录。

| 文件 | 作用 |
|---|---|
| `source/s123_four_rows.py` | 利用第四行差值与二元预算的整数入口 |
| `source/verify_s123_four_rows.py` | 完整整基直接扫 \(c\) 的整数平移，再取公共 \(d\) |
| `frozen/four-row-entry.json` | 主枚举逐阶段计数 |
| `frozen/independent-entry.rows.jsonl` | 1,199个公共平移及 \(C_4=2\) 见证 |
| `frozen/independent-entry.json` | 独立回执和哈希 |
| `manifest.json` | 原文件与冻结账本散列 |

每条独立见证的格式是 `[[D,Tr(a),r,Tr(b),s,Tr(c),t,d,g],C4,A4,B4]`，其中 \(r,s,t\) 是完整整基 \(1,\omega\) 的非有理坐标，\(A_4,B_4\) 按 \([m,n]\) 表示 \(m+n\omega\)。

定理同时排除反转支持 \(\{5,6,7\}\)，任意有限次数的整数余因子均覆盖。一般均匀和猜想仍开放。
