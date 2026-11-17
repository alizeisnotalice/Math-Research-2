# 第46轮：实二次六次核心的复现包

此包支撑 `proofs/046-quadratic-core-degree-six.md` 的有限整数证书。证明的结论是：在第43轮的实二次原始块归约中，非有理原始因子 `P` 不可能有次数至多六；它**不**解决原均匀和猜想。

包内没有浮点判定。有限搜索是证明中已经给出的严格参数界之后的整数枚举；“无余因子次数界”来自保存每一层完整状态的空前缀树，而不是深度截断的阴性搜索。

## 文件布局

- `frozen/`：取得时的程序原字节；其 SHA-256 在 `SHA256SUMS.json` 中，未为可移植性改写。
- `data/single/`：位置 `{1}` 的455个二次数域核心和终层排除账本。
- `data/linked/`：位置 `{2,3}` 的35,488母行、14个精确筛选层、两棵 HGQ 前缀树与摘要。
- `data/parent/`：独立重放得到的 14 层一字节失败位置账本、直接 AB 树、两个跨实现比较证书。
- `tools/run_046_core6.py`：不改冻结源的可移植包装器。

冻结字节的来源映射为：`frozen/parent/` 对应独立重放器与树比较器，`frozen/linked/` 对应关联 `{2,3}` 核心筛选及 HGQ 树，`frozen/single/` 对应 `{1}` 核心枚举。对应输入/结果分别放在同名 `data/` 子目录；`tools/run_046_core6.py` 只在调用者提供的空工作目录中恢复这些历史相邻路径。

`single-terminal-witnesses` 原始 gzip 为 11,164,894 bytes，超过本仓库每文件 10 MiB 的归档限制。因此它被按 JSONL 行边界分成八个 gzip 文件。`data/single/single-terminal-witnesses-split.json` 给出原始压缩文件 SHA-256 `efe650350ae2a64e597f396ed93ca5db16cbb443b83a4a4093baab63c16a1000`、解压后总行数 1,061,885、解压流 SHA-256 及每片哈希。按该清单顺序拼接**解压后的字节流**即为原日志内容；分片的 gzip 字节本身不声称与原 gzip 字节相同。

## 快速独立检查

要求 Python 3（仅标准库）。在本目录运行：

```sh
python3 tools/run_046_core6.py validate
python3 tools/verify_split_jsonl_gzip.py
python3 tools/run_046_core6.py compare-trees --out /Volumes/data/CodexTemp/T2/round046-tree-compare.json
```

第二条检查分片的解压流和行数；第三条只逐项读取已有两套树，比较 7,417 个核心、55,149 层所有接受状态与 `C` 二元性拒绝数；不生成新数学搜索。原树比较用时约 0.085 秒，实际时间依机器而异。

若要重跑一个已冻结的较重验证，而不污染本包，请显式提供空工作目录，例如：

```sh
python3 tools/run_046_core6.py replay linked-b --b 1 --workdir /Volumes/data/CodexTemp/T2/round046-replay-b1
python3 tools/run_046_core6.py replay single --workdir /Volumes/data/CodexTemp/T2/round046-replay-single
```

包装器将原字节复制到临时 `round046-parent/`，并以只读符号链接提供原脚本所期待的 `agents/...` 相邻路径；结果只写入所给工作目录。`direct-prefix` 会重建全部 7,417 棵直接 AB 树，时间和输出量更大。

## 覆盖与证据链

1. `data/single/core6_single_results.json` 和独立的 `data/parent/single-core-summary.json` 给出 `{1}` 形的 `455 / 17,372 / 638,243 / 0` 逐阶段计数。高端强制 `G4=Q4=0`，所以零个四次前缀已覆盖任意有限余因子次数。
2. `data/linked/core6_linked_conditional_ledger_corrected.json` 是 35,488 个 `{2,3}` 母行；它的文件内 ``conditional`` 标签只表示该脚本把根界作为输入，而该根界已在 proof 第3、5节解析证明。`data/linked/INVALIDATED_FIRST_RUN.md` 记录早期负数 `ceil` 漏项及其撤回。公开结论只依赖 corrected ledger、父重放与后续树，不依赖失效初版或历史 conditional ledger。
3. 14 个 `data/parent/linked-b-*-failure-codes.bin.gz` 和对应 JSON 摘要记录 22,889,760 个 `(d,e)` 项的首失败位置。`data/parent/linked-complete-crosscheck.json` 绑定与原筛选器完全一致的根筛结果：8,635,471 个正轴通过、7,417 个严格 Schur 留存。
4. `data/parent/linked-direct-trees.jsonl.gz` 与 `data/linked/prefix_tree_{0_3708,3708_7417}.jsonl.gz` 分别是直接 AB 和独立 HGQ 枚举。`linked-tree-comparison-reproducible.json` 是已运行的比较结果；上面的 `compare-trees` 命令可独立重做该轻量比较。

## 证明绑定

组包时读取的终稿为 `proofs/046-quadratic-core-degree-six.md`，SHA-256：

`10a5226a3e1e8a916d154bcfbec8e195475ffa61c669906135621b7fb2ab5831`。

父线程将本 staging 目录整体复制至 `computation/results/046-core6/` 后，必须再记录复制后的文件哈希，并重新绑定当时最终 proof 的 SHA-256；本 README 的该值只说明本次审查所读取版本，不替代最终发布绑定。
