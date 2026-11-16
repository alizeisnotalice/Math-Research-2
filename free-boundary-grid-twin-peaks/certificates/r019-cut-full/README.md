# 原条件竖切矩阵：完整独立集

对应 ../../proofs/centered-cut-covariance.md。原始整数L/H与最终协方差分子均保留；n3还由9!全标号直接枚举校准。

```sh
clang++ -O3 -std=c++17 cut_covariance_ie.cpp -o /Volumes/data/CodexTemp/T2/a022-cut-ie
/Volumes/data/CodexTemp/T2/a022-cut-ie 6
python3 verify_n3_by_permutations.py
python3 verify_cut_covariance.py
```

后一个脚本读取本目录、相邻r019-association及r018-root-moments数据，逐整数/有理数核验。程序写小摘要，重新生成全枚举时请复制脚本到临时目录或保留当前证书。

复现时请先确认 data 挂载，将 TMPDIR/TMP/TEMP 及编译输出设到 `/Volumes/data/CodexTemp/T2`；Swift/Clang 模块缓存为 `/Volumes/data/CodexTemp/cache`。环境：Python 3.15，Apple Clang 21，Swift 6.4.2/Metal。源码相对路径调整不改变数学核；原运行源定位保留在本轮报告，归档哈希见 manifest。有限核验不外推到所有 n。
