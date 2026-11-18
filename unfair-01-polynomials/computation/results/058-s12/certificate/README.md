# 第58轮位置 \(\{1,2\}\) 排除证书

[中文定理和覆盖证明](theorem-zh.md)。

本包无需下载数据或安装库。Python 3.13.18 标准库重算112个实二次整环域，两种实现分别用半整数根式和完整整环基。第二套直接遍历允许的 \(d,f,g\)，不借用第一套的第五行差值定位和通常分支中的 \(f\) 定位。运行时会验证源文件及冻结基准的散列，然后重新生成全部账本、逐母式比较第一拒绝原因，并与冻结结果逐字节对照。

```sh
/opt/homebrew/bin/python3.13 /path/to/package/replay.py \
  --out /Volumes/data/CodexTemp/T2/round058-s12-replay --workers 4
```

选择包外不存在的输出目录。可把解释器换成可用的 Python 3，必须启用断言。每次输出含两个程序的 stdout 和 stderr、三份新账本、所有子进程退出回执与文件散列。

| 文件 | 含义 |
|---|---|
| `source/s12_exact_entry.py` | 根式坐标的完整枚举及严格根筛 |
| `source/verify_s12_integral_basis.py` | 整环坐标的独立枚举及首拒见证 |
| `frozen/full-entry.json` | 域、行数、反转轨道账本 |
| `frozen/full-entry.rows.jsonl` | 全17项及首拒码 |
| `frozen/independent-entry.json` | 正点评估见证与范数 Schur 首拒次数 |
| `manifest.json` | 可移植源及输入文件 SHA-256 |

每个母式记录为 `[D,Tr(a),r,Tr(b),s,c,d,e,f,g]`，其中 \(\omega=(1+\sqrt D)/2\) 当 \(D\equiv1\pmod4\)，其余为 \(\sqrt D\)。\(r,s\) 为两个 \(\omega\) 坐标。独立账本的点评估 `value_basis=[u,v]` 表示清分母后的 \(u+v\omega\)，`side=0` 表示 \(P\)，`side=1` 表示 \(\bar P\)。外环见证包含整数范数升幂系数 `H` 和首失败次数。

本包只排除八次实二次原始母式的两个反转支持 \(\{1,2\}\) 与 \(\{6,7\}\)。一般猜想和其余支持仍开放。
