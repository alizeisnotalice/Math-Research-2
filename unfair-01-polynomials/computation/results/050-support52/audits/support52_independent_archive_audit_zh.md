# total 52、cap 6 支撑归档的独立审查

审查日期：2026-11-17。审查者只读了公开归档及其冻结输入；未启动远端、GPU 或生产任务。

## 审查对象与结论

审查对象为

`/Volumes/data/CodexTemp/T2/unfair-01-polynomials/sync/unfair-01-polynomials/computation/results/050-support52/`。

**接受如下精确而有限的结论：** 对总次数 (52)、两端点归一化、共同指数公因子为 (1)、每个支撑至多六点的支撑对，在通常的交换及同时反转约化后，归档完整枚举了每一个代表；其冻结的 singleton-unit 证书分类没有留下 `unknown`。这是一项有限支撑排除/公平分类结果，不能推出任意次数、任意支撑大小或原始一般猜想。

这不是“所有总次数 52 支撑”不加限定的结论：本次只枚举 cap (=6) 且 primitive（共同指数公因子 (1)）的对。非 primitive 对可按共同指数约化到较小总次数，但是否已由别轮覆盖不是本归档单独证明的事项。

## 实际复跑

运行：

```sh
python3 verify_support52_archive.py \
  --package archive \
  --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round050-support52-audit/archive-reverify.json
```

实际输出为 `passed_parent_full_archive_hash_binding_and_support_coverage`，耗时 6.358 秒；输出 SHA-256：

```
95afcc02dc9cebb6e5a4093b25ffbacfaefc2313058a335dacd244680620cdc2
```

该次复跑重新生成 (83,681) 个 (I)-shard 的索引，逐个核对缓存回执、事件哈希、输入/输出压缩哈希、绑定及总计数。聚合数为：

| 项目 | 数量 |
| --- | ---: |
| raw | 2,490,683,871 |
| canonical | 1,185,395,011 |
| excluded | 1,178,648,372 |
| fair | 6,746,639 |
| unknown | 0 |

本次所绑定的 manifest SHA-256 为 `8087c3ff5edea0c65edb3849ab0371cac014c4313d6f3d00d519d6b47d56ff51`；独立 chunk 语义验证器 PIN 为 `a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b`，并已由共享归档中的压缩源重新验 hash。

## 覆盖与半平衡 (26+26) 的核查

固定 (m\leq26)，端点归一化、至多六点的 (I\subseteq\{0,\ldots,m\}) 的数目为

\[
 \sum_{k=0}^{\min(4,m-1)}\binom{m-1}{k}.
\]

因此总 shard 数为

\[
\sum_{m=1}^{26}\sum_{k=0}^{\min(4,m-1)}\binom{m-1}{k}
=\sum_{k=0}^{4}\binom{26}{k+1}
=26+325+2600+14950+65780=83681,
\]

与实际事件、缓存记录及验证器断言相同。对于 (m<26)，交换因子把支撑对送到未枚举的 (m>26) 一侧，故只取 (m\leq26) 正好各取一次。对于半平衡情形 (m=n=26)，代码先按 `(len(I),I)` 与 `(len(J),J)` 排序以商去交换；再把同时反转后的两支撑同样排序，保留原对当且仅当其排序对不大于反转排序对。因同时反转是二阶对合，这恰在每个“交换 + 同时反转”轨道取一个代表（若相等则该代表固定），没有遗漏或重复。

`gcd(*I,*J)!=1` 的丢弃是共同指数缩放约化，而非静默的随机过滤；此项须与主稿中已证明的共同指数公因子约化一同引用。

## 语义、原始数据和可复现性边界

冻结的 `audited_chunk_verifier.py` 并不导入搜索器：它独立生成 (J)、独立施行交换/反转/公因子约化，重算每个对的 singleton-unit 分类，且逐行比较原始 chunk。父运行器已对 83,681 个 total-52 chunk 调用该验证器；每个归档回执绑定事件 SHA、计划、冻结生成器、payload、启动清单、终端回执及上述 verifier PIN。`raw52.json` 还报告对 4,528,314,671 bytes 压缩原件（65,378,298,582 bytes 解压原件）的逐 chunk 哈希/行数复核，无坏块。

需要准确措辞：6.358 秒的公开归档复跑是**所有语义回执、覆盖与字节绑定的重新核验**，不是对 65GB 原始 chunk 再进行一次语义重算；语义重算是在已绑定的父运行器中完成。当前小归档未包含那 65GB 原始 chunk，只包含其 `raw52` 审计摘要。因此公开包应同时发布/索引 `raw52_audit.py`、`raw52.json`、生成器/runner 及其哈希，并明确原始 chunk 的冻结位置和哈希清单；这会使“可重建而非把成功标记当证据”的链条可审查。此为文档/索引补项，不否定已通过的有限覆盖及现有哈希链。

## 审查决定

无发现会导致 total-52、cap-6、primitive 有限覆盖或 (26+26) 代表约化失效的漏洞。可在完成上述原始数据索引披露后，将此归档作为该**限定范围**的计算证据；不得据此登记一般猜想已解决。
