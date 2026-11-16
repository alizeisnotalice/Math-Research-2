# r018核的独立高粒子n8对照

P16384,R32,L16,8线程，基础seed0xA0222100，2016条完整重复。源是实际执行的thread_local版本，未加入后来归档入口guard；本次参数均在其安全范围、TRACE未启用。该差异不被隐瞒为同一文件哈希。

```sh
clang++ -O3 -std=c++17 -pthread b1.cpp -o /Volumes/data/CodexTemp/T2/a022-high-p
/Volumes/data/CodexTemp/T2/a022-high-p 8 16384 32 16 8 /Volumes/data/CodexTemp/T2/a022-high-p.csv 0xA0222100
python3 ../r018-cpu/audit_cpu_csv.py n8-p16384-L16.csv
```

root-summary用各秩均值之和取比，分层delta SE不是重复比值均值的SE。所有数值均仅诊断；本批已完成退出。

复现时请先确认 data 挂载，将 TMPDIR/TMP/TEMP 及编译输出设到 `/Volumes/data/CodexTemp/T2`；Swift/Clang 模块缓存为 `/Volumes/data/CodexTemp/cache`。环境：Python 3.15，Apple Clang 21，Swift 6.4.2/Metal。源码相对路径调整不改变数学核；原运行源定位保留在本轮报告，归档哈希见 manifest。有限核验不外推到所有 n。
