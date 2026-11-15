# 计算复现

所有脚本使用Python标准库；不需要安装numpy/sympy。精确计数使用int/Fraction或C++ unsigned __int128；n≤5限制不可擅自放大。随机搜索不承担渐近证明或混合认证。

在 `/Volumes/data/CodexProjects/T2` 工作。依赖data写入前确认其已挂载。临时输出根为 `/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks`，C++共享缓存为 `/Volumes/data/CodexTemp/cache/clang`。

## 核验已归档证书

```sh
python3 -B free-boundary-grid-twin-peaks/computation/verify_results.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_insertion_n3.py
python3 -B free-boundary-grid-twin-peaks/computation/compute_global_geodesic_p.py
python3 -B free-boundary-grid-twin-peaks/computation/thin_cycle.py
```

## 重建精确计数

```sh
python3 -B free-boundary-grid-twin-peaks/computation/exact.py --max-n 4 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/exact-rerun.json
```

每个C++源用下列方式编译（把moment_dp换为cut_dp、pair_dp、geodesic_dp、merge_tail或conditional_walk）：

```sh
TMPDIR=/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -std=c++17 -O3 -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang free-boundary-grid-twin-peaks/computation/moment_dp.cpp -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/moment_dp
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/moment_dp 5
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/geodesic_dp 5
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/merge_tail 5 8
```

merge_tail的第二参数t给严格尾事件M>t；将t遍历2至25即可复现完整尾部。峰对与切割扫描：

```sh
python3 -B free-boundary-grid-twin-peaks/computation/pair_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/pair_dp --workers 3 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/pairs-rerun.json
python3 -B free-boundary-grid-twin-peaks/computation/cut_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/cut_dp --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/cuts-rerun.json
```

## 随机实验

```sh
python3 -B free-boundary-grid-twin-peaks/computation/search.py --steps 200000 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/search-rerun.json
python3 -B free-boundary-grid-twin-peaks/computation/walk_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/conditional_walk --steps 5000000 --workers 3 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/walks-rerun.json
python3 -B free-boundary-grid-twin-peaks/computation/walk_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/conditional_walk --steps 5000000 --workers 3 --reflect-every 10 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/walks-reflect-rerun.json
python3 -B free-boundary-grid-twin-peaks/computation/walk_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/conditional_walk --steps 20000000 --workers 3 --reflect-every 10 --top-block-every 4 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/walks-block-rerun.json
```

种子由源码固定并逐链记录。walk前20%作烧入，之后每n²步记录一次；保留前后半段直方图和独立核验的最佳候选。提议包括固定分布的标签/位置交换、可选条带对合，以及可选含最高峰2×2块的均匀置换。块由最高峰位置选择时，接受比为该点原有可选块数/新位置可选块数；其他提议对称。所有更新均要求最终恰有两峰。均匀条件律平稳不等于不可约或混合已证。

浮点仅用于退火接受率与探索摘要；精确证书不依赖这些值。不同C++标准库实现的shuffle可能改变随机轨迹，实际编译器版本已记录。旧阶段输出保留当时二进制哈希；最终源码通过相应参数保留三阶段算法，性能计时不是复现相等项。

## 第2轮新增核验

`weighted_promotion.py` 完整枚举3×3，验证所有提升关系元、逐输出归一化、精确逆像及有理对数证书。`promotion_dp.cpp` 接收边长n（2≤n≤5），以第一触邻递推给出总提升质量的精确分子/分母；n5内存1 GiB。编译缓存沿用上节，命令和实际结果在rounds/002.md与manifest-r002.json。

## 第3轮新增核验

```sh
python3 -B free-boundary-grid-twin-peaks/computation/verify_r003.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_shortest_lift.py
python3 -B free-boundary-grid-twin-peaks/computation/beta_fiber_n3.py
python3 -B free-boundary-grid-twin-peaks/computation/geodesic_fiber.py --n 4 5 6 7 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/geodesic-fiber-rerun.json
```

按上节同样编译singleton_dp.cpp与beta_walk.cpp；前者n5内存1 GiB，命令为`singleton_dp 5`。后者新增任一峰2×2块和β观测，复现命令：

```sh
python3 -B free-boundary-grid-twin-peaks/computation/walk_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/beta_walk --steps 20000000 --workers 3 --reflect-every 10 --top-block-every 4 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/beta-walks-rerun.json
```

此处沿用驱动参数名top-block-every，但beta_walk中从两峰等概率选择。给定块B，选择概率正比于Σ_{峰v∈B}1/b(v)，b(v)是含v的2×2块数；接受比为输出/输入的此项之比。仅理想随机核的平稳性有解析依据，实际浮点/伪随机实现及有限轨迹不提供认证概率。全部最佳标号以独立峰判定器再次核验。
