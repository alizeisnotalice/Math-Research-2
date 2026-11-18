# 第52轮复现

运行环境为 ARM64 macOS、Python 3.15.0、Apple clang 21.0.0 和现有 Swift Metal 后端。
在 data 卷选择空目录保存新产物。下列命令从本证书目录运行，SCRATCH 指向空目录。

## 解析核验与全图映射

```sh
python3 -B analysis/audit_group_max.py | cmp - analysis/audit_group_max.stdout
python3 -B analysis/audit_exclusive.py | cmp - analysis/audit_exclusive.stdout
python3 -B roots/audit_roots.py | cmp - roots/audit_roots.stdout
python3 -B construction/promote_seams.py "$SCRATCH/seam-gate.json"
cmp "$SCRATCH/seam-gate.json" construction/map-gate.json
```

核验项目包括真实小规模后缀的精确整数、固定坐标几何、阈值积分、全 \(9!\) 输入的
指定对角角点映射，以及带固定种子的额外正逆重放。后两种构造样本只核验映射合法性。

## 十二阶证书

```sh
python3 -B construction/extend_n12.py "$SCRATCH/numerator.json"
cmp "$SCRATCH/numerator.json" construction/result.json
python3 -B t5/refine_t.py ../r051/n12/denominator-final.json \
  ../r051/n12/n12-conditional-final.json "$SCRATCH/denominator.json"
cmp "$SCRATCH/denominator.json" t5/r052-denominator-refined.json
python3 -B t5/check_transpose.py
python3 -B t5/join_bounds.py | cmp - t5/combined.json
```

如需重算原始模数通道，将 t5 中的 rect_protected.cpp、run_t5.py、pratt-61bit.json
复制至计算目录。设共享编译缓存为 /Volumes/data/CodexTemp/cache/clang，以
`clang++ -std=c++20 -O3 -Wall -Wextra -Werror` 编译为 rect_protected。
本机运行 `python3 -B run_t5.py local 3 4 --bottom`，MacBook 运行
`python3 -B run_t5.py macbook 0 1 2`。两个批次各自读取本机编译产物，内存独立。
再把两个输出目录放在同一个计算目录，使用 refine_t.py 的 `--shards` 指定它。
对归档复算时直接使用现有模数通道、哈希和退出回执。测时、PID 与编译产物哈希会变化。
上区域另本机重跑转置 `./rect_protected 5 12 R 0 --max-table-mib 3000 --progress`，
三张模数直方图和全部层大小与 MacBook 原方向逐项相同。

## 本机十六阶 Metal

长短两批均使用现有第26轮 Swift 源，\(n=16\)、131,072 粒子、12 重复、\(\tau=2\)、
种子域 9e52128000000001，分层重采样。移动步数分别为128和16。
从同一冻结源码用 `swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift`
编译至独立目录，在长短两目录分别运行 local-metal/run.py 和 local-metal/control/run.py
的副本。验证器检查全部种子、原始文件和源码哈希、退出回执及整数矩关系，
并以质量加权比较两批的有限诊断。较早第28轮的账本保留为额外参照。

```sh
python3 -B local-metal/verify.py | cmp - local-metal/summary.json
```
