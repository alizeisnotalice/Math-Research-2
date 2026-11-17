# Independent combined total-46 verification

This archive verifies exactly the endpoint-normalized, at-most-six-term
support search slice at total degree 46.  It combines the reboot-interrupted
old producer's 1,998 verified chunks with the new missing-only resume's 42,553
chunks.  It does **not** claim a proof for arbitrary degrees or supports.

The independently enumerated union has 44,551 support-index chunks, with no
overlap or omission.  The independently recomputed totals are:

```text
raw       787855938
canonical 372691286
excluded  370540360
fair        2150926
unknown           0
```

`combined46-full-report.json.gz` holds the complete semantic report; its gzip
SHA-256 is `e502656f6707a1c9b916c79c9ba03e511342374231a332eba0551ba69bdb9ed6`
and its uncompressed SHA-256 is
`3563408d7502205876aea4de220cfed510ce0fe4f63fcf1b08c049af7d761748`.
`combined46-chunk-bindings.json.gz` is the per-chunk hash/count index
(gzip SHA-256 `016037a8d95aaf36cf58ccaa508c10fea8bd1f3901af3978848ff3bc0b586f0f`).

The inputs are bound by the old event SHA-256
`10ec3f5465b3a3ab8d6d5d73922967525267905d7770a4eca3425bab4f8bec78`, new
event SHA-256 `869499025260a928fff2b1a4f955e12d4517c9b209d03185265d50b4bc39e587`,
frozen classifier SHA-256
`d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710`,
and independently audited verifier SHA-256
`a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b`.

The negative tests rejected a deleted chunk, a duplicated chunk, and a
modified singleton-unit witness even after recomputing that fixture's gzip and
raw hashes.  Source code copies are under `sources/`.

## 本轮主线程复核与范围

独立逐项语义核验实际用4工作进程耗时1925.6396秒，全部44551片通过；主线程另用 check_combined46_coverage.py 全解压20,266,444,849字节账本并重算覆盖、来源绑定与哈希，47.9389秒通过。后者不是再次进行全部数学语义判定；语义结果来自前者的完整重算。两者共同支持新增总次46层，与既有≤45合并得每侧≤6点、总次数≤46严格公平。

sources 中保留独立验证器与冻结生产源；new-run 和第39轮039-support45的旧流保存重建入口。完整可再生大账本留在 CodexTemp 的原始 run 目录，压缩报告不代替逐项原输入。需要复现全部运行时，按各脚本 --help 指定 old/new run、plan 与输出路径，先由归档生产源重建分片，再运行独立验证器。parent-archive-binding.json 说明原始 pretty JSON 与归档 canonical JSON 为同一对象，二者字节哈希不可混用。

本机完整语义复现命令（输入大账本仍需保留；下列只写独立新输出）：

```sh
/Users/zzh/.venvs/main/bin/python sources/verify_combined46.py \
  --old-run /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-support46/run \
  --new-run /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round040-macbook-resume46/run \
  --plan /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/resume46_macbook40/plan/resume46-plan.json \
  --frozen sources/frozen_support.py --remote-resume sources/remote_resume46.py \
  --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/recheck46/report.json \
  --cache /Volumes/data/CodexTemp/T2/unfair-01-polynomials/recheck46/cache --workers 4
```
