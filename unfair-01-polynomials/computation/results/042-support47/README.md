# Independent total-47 certificate

This archive certifies only the fixed total-degree 47, cap-6 support layer:
44,551 chunks independently enumerated, hash checked, and semantically
recomputed.  It passed with `unknown=0`; this is not a general theorem.

`total47-full-report.json.gz` and `total47-chunk-bindings.json.gz` are each
about 2 MiB.  The report SHA-256 is
`316aaca8465d79f124b78c26f48cd16d3a9d6e70dbfa35ed732725465ed7464d`.
The source event ledger is
`/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round040-macbook-support48/run/events.jsonl`;
plan, frozen classifier, and producer are respectively
`/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/support48_macbook40/plan/support48-plan.json`,
`.../plan/frozen_support.py`, and `.../remote_support48.py`.

Complete replay command:

```sh
/Users/zzh/.venvs/main/bin/python sources/verify_support48.py 47 \
  --run /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round040-macbook-support48/run \
  --plan /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/support48_macbook40/plan/support48-plan.json \
  --frozen /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/support48_macbook40/plan/frozen_support.py \
  --remote /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/support48_macbook40/remote_support48.py \
  --out replay47.json --cache replay47-cache --workers 4
```

The same production payload also contains total 48, but its independent
verification is still running and is intentionally absent from this archive.

## 主线程覆盖与归档绑定

独立语义核验已检查44551片，2399.311秒，未知0；其报告与搜索器分离。主线程另运行 check_total47_coverage.py，用独立组合枚举重建全部I输入，重算每片压缩及未压缩hash、行数、原始J计数与报告绑定：完整解压24,932,589,389字节，457329424行，55.894秒退出0。该步骤核内容与覆盖，不冒充再次语义重算。

production中保存完整47/48生产事件和payload压缩副本；本归档只认证47。总次48仍在独立验证，不纳入严格界。完整可再生原始分片位于README上列run/t47，可由冻结生产源与计划重建。压缩报告采用canonical JSON；parent-archive-binding明确与原pretty JSON完全相同对象，按indent=2,sort_keys=True加末尾换行可还原原始hash。
