# r015 资源与检查点

- 机器：macOS 27.2，物理内存 17179869184 bytes，`hw.ncpu=12`。
- 编译器：Apple clang 21.0.0；GMP 来自 `/opt/homebrew/opt/gmp`。
- 线程：完整 n=5 用 4 个 CPU 线程；原因是五个全图 `unsigned __int128` 数组已经
  占 (5\cdot2^{25}\cdot16=2.50\) GiB，且同机还有 GPU 驱动和其他用户进程。
- 检查点：`n5-benchmark-1.json`（一项、(2^{24}) 压缩状态）和
  `n5-benchmark-100w4.json`（100 项）是全量前的真实运行输出。
- 完整 n=5 的最终核验运行：24.67 wall seconds，69.88 CPU seconds，maximum RSS
  4364959744 bytes，peak footprint 5100146672 bytes，`swaps=0`。程序进程已退出。
