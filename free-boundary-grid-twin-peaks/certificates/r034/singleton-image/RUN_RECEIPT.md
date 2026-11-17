# R037 独立枚举运行收据

- 日期：2026-11-17（Asia/Shanghai）。
- 解释器：`Python 3.15.0`。
- 实际命令：

  ```sh
  time nice -n 10 python3 \
    audit_singleton_image_n3.py
  ```

- 退出码：0；实际计时：`1.95s user 0.01s system 46% cpu 4.220 total`。
- 输出：`audit-singleton-image-n3.json`，其 SHA-256 为
  `d98278424406c968c017da1d124f51be3468e83a600eb2011697ea4f062737f5`。
  它实际构造所有一峰源及 eligible $p$ 的提升、以 `Counter` 数每个输出的逆像、逐
  输出核对这个实际数等于 $R=b-l-1$，再独立扫描全体 top-two 词并按第一次合并判定
  low-singleton。没有把先前 R 公式的调和求和当作像枚举的替代品。
- 源码 SHA-256：
  `f15970d47a8ed7621d0199ceced199c83695f52f5907675aa77bce14d1744196`。
- 运行结束后没有该 Python 进程；未用 SSH、GPU 或大内存计算。
