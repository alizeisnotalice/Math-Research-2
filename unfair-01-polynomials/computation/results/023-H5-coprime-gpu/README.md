# 第23轮 H 五项母式互素替换：最小 GPU 归档包

本包保存 `H=1+x+x^2+x^5+x^6` 的 11 组数值根组筛选账本，供复现和
独立审核。它不提供一般猜想的反例或严格有限排除；GPU 阴性结果只用于数值
探路。相应 H 双替换族的公平性必须由单独的 carry 解析证明承担。

本临时归档的最终公开位置由 `archive-map.json` 固定为
`unfair-01-polynomials/computation/results/023-H5-coprime-gpu/`；提交时应将
整个包按该子目录原样复制，而不是保留对 `agents/` 路径的依赖。

## 内容和大小边界

`checkpoints/` 是 11 份原始持久化账本；每份均小于 10 MiB，未经截断。
`inputs/` 是逐组根数据和精确二元目标；内容与原输入逐字节相同。`receipts/`
含搜索后独立 MP80 审计回执。`scripts/stream_coprime_ds_gpu_v1.py` 保持原字节，
SHA-256 为
`4ecb8f50db982f90a2467a5cd25fbe109d6c98d11fa26439a2d4dfed46bcafee`。

选项缓存没有归档：它们是可再生的数值根组表，且个别文件很大。`
scripts/prepare_coprime_inputs.py` 和 `scripts/verify_coprime_layout.py` 可用
MP100 重建和验证它们。该准备器固定为本 H 五项母式、三个母根块、
`1<a<b` 且互素的输入，不适用于本轮以外母式。

## 重建缓存

以下命令以完整 checkpoint 的 `(a,b)=(5,8)` 为例；临时项须写在 data 盘。输入中的 SHA-256
必须与 `archive-map.json` 相符，才能用于同一 checkpoint。

```sh
PY=/Users/zzh/.venvs/main/bin/python
ARC=unfair-01-polynomials/computation/results/023-H5-coprime-gpu
CACHE=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/rebuild-cache/a5b8-options.json
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  "$PY" "$ARC/scripts/prepare_coprime_inputs.py" --a 5 --b 8 \
  --out-input /Volumes/data/CodexTemp/T2/unfair-01-polynomials/rebuild-cache/a5b8-input.json \
  --out-options "$CACHE"
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  "$PY" "$ARC/scripts/verify_coprime_layout.py" \
  --input "$ARC/inputs/a5b8-input.json" --options "$CACHE" \
  --output /Volumes/data/CodexTemp/T2/unfair-01-polynomials/rebuild-cache/a5b8-layout.json
```

The regenerated input should be checked byte-for-byte or by the recorded
input hash before an archived checkpoint is paired with it.  Numerical root
preparation is not a root-isolation certificate.

## Search and independent checks

The CPU-only layout/resume test is:

```sh
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  "$PY" "$ARC/scripts/test_coprime_layout_resume.py"
```

A new MPS screen is invoked explicitly, for example:

```sh
"$PY" "$ARC/scripts/stream_coprime_ds_gpu_v1.py" \
  --options "$CACHE" --roots-input "$ARC/inputs/a5b8-input.json" \
  --checkpoint /Volumes/data/CodexTemp/T2/unfair-01-polynomials/rebuild-cache/a5b8-checkpoint.json \
  --batch-size 131072 --checkpoint-seconds 30 --checkpoint-interval-batches 256
```

For a complete checkpoint, the frozen independent verifier must be used
without `--allow-partial`:

```sh
"$PY" "$ARC/scripts/verify_coprime_result.py" \
  --checkpoint "$ARC/checkpoints/a5b8.json" --options "$CACHE" \
  --roots-input "$ARC/inputs/a5b8-input.json" \
  --engine-source "$ARC/scripts/stream_coprime_ds_gpu_v1.py" \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/recheck-a5b8.json \
  --mp-dps 80 --tolerance 1e-7
```

`a5b9` and `a3b10` are intentional partial prefixes after analytic closure.
Their final 64 saved rows lack the searcher's CPU float64 copy, so the frozen
verifier correctly refuses them. The explicit-path supplement
`audit_missing_saved_rows.py` re-decodes every retained row and independently
reconstructs it at MP80; it does not change a checkpoint or call the prefix
complete. Both archived-copy rechecks were actually run, with receipts
`receipts/a5b9-archive-mp80-recheck.json` and
`receipts/a3b10-archive-mp80-recheck.json`.

## Exit and proof boundary

`exits/H5-compute-final.json` is the final consolidated process-reaping
receipt. It records which queues completed and which were sent SIGINT after a
separate proof closed their family. It must not be read as a full enumeration
claim for those partial prefixes. `archive-map.json` specifies the intended
topic-relative destinations and each file's original source.
