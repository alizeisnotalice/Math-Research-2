# 实二次五次原始块完整排除证书

完整解析覆盖见 `../../../proofs/045-quadratic-core-degree-five.md`。不是给G/Q设次数49上限：由根界得有限母式全集，404必要母式的全部自由前缀树在49行以内死亡，任意有限G/Q补零后必属于树，因此排除任意次数。

## 实际核验

- 搜索：3,981,393整数格；判别式/嵌入后584,820；七正有理点后176,049；严格双单位圆递推后404。
- 主线程 `verify_core5.py` 独立重建H=F²+Tx²F+Nx⁴，以降幂整数递推检查每项首个失败见证、完整覆盖及404列表，实际8.583秒通过。
- 主线程 `direct_prefix.py` 直接在Z[a]重算A*B全部卷积，不导入搜索器，404树全部死亡，最晚49层、最大7状态；实际0.028秒。
- 独立 `hgq_prefix.py` 用H(GQ)重算，明确含h0 Gk Q0/h0 G0 Qk以及Pk G0端项。所有404逐层计数相同。父再核完整3,107层所有接受状态及全部终层拒绝项一致，回执 `cross-implementation-state-check.json`。
- Python实际3.13.18；所有判定为整数计算。深度64与状态100000是程序防护阈值，本次404均提前为空，没有未知、超时或截断项。

## 文件和复现

`enumerate_core5_candidates.py` 是筛选器，`core5_summary.json`及两份压缩JSONL包含完整覆盖和排除见证；解压后文件hash须与summary一致。`core5_survivors.json`为全部404输入。

`direct-prefix-trees.json.gz`保存父每层全部父前缀、允许延伸和拒绝卷积系数。`hgq-all-layer-certificates.json.gz`把第二实现各core/layer与death JSON的语义内容合并为按相对路径索引的无损字典；原逐gzip的hash仍在`hgq-summary.json.gz`，合并包不冒充原gzip字节。最终证明审查绑定见`final-proof-audit.md`。

本次原始目录为 `/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/coprime_stretches/round045-core5`，父实现目录为 `/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round045-parent-core5`。使用 `/Users/zzh/.venvs/main/bin/python` 依次运行 enumerator、父verify_core5、父direct_prefix和prefix/hgq_prefix即可复算。公开包异地复现时，先将两JSONL解压并把父脚本ROOT/SOURCE常量指向解压目录；将hgq_prefix放在该目录的prefix子目录。输出另放副本，避免覆盖已冻结证书。无需GPU或数值库。
