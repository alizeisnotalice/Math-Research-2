# 第43轮复现

将 PKG 设为本目录，OUT 设为数据盘上的全新目录。Python 3.15 标准库与 Apple clang 21 已实际运行。

```sh
env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
  python3 -B "$PKG/replay.py" "$OUT"
```

入口读取第42轮区域计数门槛和已归档的精确七阶双峰分母。重建32项小门槛、4项七阶代表与切割递推核验。完整128条七阶轮廓在 MacBook 上运行。每条轮廓的两侧计数、交错数及归一化都保存在公开账本。独立子集递推采用 `region_subset.cpp`，编译选项同 replay，参数为 `3434343 OUTPUT.json` 或 `4343434 OUTPUT.json`，峰值驻留内存分别约270和539 MB。

MacBook 的三组混合抽样输入和结果摘要在 mixture 子目录。用 `clang++ -std=c++20 -O3` 编译 `exact_curved_mixture_sampler.cpp`，然后运行：

```sh
python3 -B "$PKG/mixture/run_mixture_remote.py" "$PKG/mixture/n7.json" \
  "$BINARY" "$NEW_DIR" --samples 65536 --master-seed 2026111803 \
  --workers 8 --cache-mib 16 --chunk-size 256
```

n3与n5的样本数为8192和16384，种子为2026111801和2026111802。MacBook 已收齐全部作业和退出回执。抽样原词留在项目本地审计材料，公开归档保留输入、版本和摘要。将 `verify_r077.py` 的第一个参数设为七阶抽样输出目录（含 raw 和 summary.json）即可复核全部重数直方图。实际调用记录见 `r077-samples-audit.json`。有限 PRNG 输出用于估算所列事件并集。
