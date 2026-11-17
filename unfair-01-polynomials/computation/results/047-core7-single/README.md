# 第47轮：七次单位置 `{1}/{6}` 排除的复现包

本包服务于实二次原始块模型
\(A=GP,B=Q\bar P\)，而非一般均匀和猜想。它严格排除七次 \(P\) 的非有理支持
`{1}` 及其同时反转 `{6}`；其余五个反转支持类仍未解决。

## 内容与来源角色

- `proofs/`：七次支持/根界归约与 `{1}` 分支完整证明。
- `frozen/parent/`：独立的 HGQ 前缀核验器及逐核心 AB/HGQ 比较器原字节。
- `frozen/independent/`：独立直接 AB 低端筛和其680核控制器原字节。
- `data/parent/`：680条 HGQ 记录、汇总、SymPy 末行化简和 AB/HGQ 全量比较结果。
- `data/independent/s1-ab-records-680.jsonl.gz`：680个直接 AB 原始记录的逐条无损合并。每一行保存原文件名、原始 UTF-8 字节和该字节的 SHA-256；没有提交680个零散文件。
- `data/audit/`：独立输入盒、账本、终端 Newton 算术审计。

所有程序只使用 Python 标准库；实际原运行以严格整数/平方根比较完成，没有浮点判定。

## 快速核验

在本目录执行：

```sh
python3 tools/validate_manifest.py
python3 tools/verify_merged_s1_records.py
mkdir /Volumes/data/CodexTemp/T2/round047-s1-records
python3 tools/extract_s1_records.py /Volumes/data/CodexTemp/T2/round047-s1-records
python3 frozen/parent/compare_s1_complete.py \
  --parent data/parent/s1-hgq-all680.jsonl \
  --records /Volumes/data/CodexTemp/T2/round047-s1-records \
  --source frozen/independent/s1_lowprefix_screen.py \
  --out /Volumes/data/CodexTemp/T2/round047-s1-comparison.json
```

最后一条逐项比较680个 \((T,N,\Delta)\)：输入索引、`c2/c3/c4/prefix4/c5/prefix6` 的所有计数和每个存活状态必须完全相同。已归档结果
`data/parent/s1-complete-comparison.json` 的状态为
`passed_all_680_input_triples_stage_counts_and_survivors`。

## 证据范围

根包络给 `T∈[-16,16]`、`N∈[-72,72]` 和严格判别式条件，恰680个域核心。系数盒为

\[
c_2,c_5\in[-29,29],\quad c_3,c_4\in[-52,52],\quad c_6\in[-8,8].
\]

高端预算强制 `G5=Q5=0`。完整筛仅留下
`(T,N,Δ)=(-2,-2,12)`、`(c2,c3,c4,c5)=(1,4,1,5)`；它被第五 Newton 和的精确不等式排除。`c6` 未列举并不产生缺口，因为矛盾已在五次前缀出现。

证明稿 SHA-256：

- `proofs/degree7_single_support_exclusion_zh.md`：`ade46954b8970be6daf41e5f9532ec9eaa80a589c1e9172002d53aba938453b3`
- `proofs/degree7_supports_root_bounds_and_rank_zh.md`：`393bdb2168bd13b34d14c1f614f29c5977603edf8cf817649365960cf2bd05e9`

本包不声称排除其他七次支持，也不替代后续精确候选搜索。
