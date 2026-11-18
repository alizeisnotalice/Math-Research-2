# 第44轮复现

将 PKG 设为本目录，OUT 为数据盘上的新目录。Apple clang 21 与 Python 3.15 已用于本机复算。MacBook 的 n8 批使用 Apple clang 21、Python 3.9.6。

```sh
mkdir "$OUT"
env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
  clang++ -std=c++20 -O3 -Wall -Wextra -Werror "$PKG/cut_region.cpp" -o "$OUT/cut_region"
export A022_R044_REGION_BINARY="$OUT/cut_region"
python3 -B "$PKG/audit_interface.py"
python3 -B "$PKG/certify_disjoint.py" --workers 6 --out "$OUT/disjoint.json"
python3 -B "$PKG/certify_parity.py" "$OUT/parity.json"
python3 -B "$PKG/run_n8.py" "$OUT/n8" --workers 3
```

前两组有限账本和接口审计已从最终公开源码重放。n8 按 MacBook 生产回执保留全部16项计数、编译版本和总和。公开脚本只改了输入与输出路径，以便不覆盖原账本。n8 单进程最大驻留约790 MB。16个计数进程均已回收。旧 n4–6 分母是包含全峰数和原距离表的原始整数输出，n7 分母沿用第5轮公开证书。
