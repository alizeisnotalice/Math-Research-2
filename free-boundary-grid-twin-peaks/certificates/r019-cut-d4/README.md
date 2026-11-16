# D4切矩阵：n7

固定竖切不具有D4不变性。每个代表显式算八变换统计和，再乘轨道大小/8；在整数系数层检查整除。n4..6与完整表全矩阵一致；n7总和与上一轮点积矩一致。

```sh
clang++ -O3 -std=c++17 -pthread cut_covariance_d4.cpp -o /Volumes/data/CodexTemp/T2/a022-cut-d4
/Volumes/data/CodexTemp/T2/a022-cut-d4 7
python3 verify_n7_d4_cut.py
```

n7运行使用8线程，78.22秒、最大RSS4,494,884,864 B、峰memory footprint6,371,038,656 B，零swap并退出。两种内存口径不同。本轮资源上限为6GiB footprint。

复现时请先确认 data 挂载，将 TMPDIR/TMP/TEMP 及编译输出设到 `/Volumes/data/CodexTemp/T2`；Swift/Clang 模块缓存为 `/Volumes/data/CodexTemp/cache`。环境：Python 3.15，Apple Clang 21，Swift 6.4.2/Metal。源码相对路径调整不改变数学核；原运行源定位保留在本轮报告，归档哈希见 manifest。有限核验不外推到所有 n。
