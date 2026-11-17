# 两个有限搜索候选

这两个保存词是正在进行的较大批次中已经完成并逐条核验的候选，不是全局最小值，也不是条件概率样本。
随机种子、10^7提案、根、初态、接受数与审计间隔1000见同名JSON和词首行。
冻结搜索源SHA与远端二进制SHA见SOURCE_PROVENANCE.json；标准mt19937_64，温度0.05。

仅检查保存词（标准库Python，有理数）：

```sh
python3 -B verify_psi.py n8-bfs-center-s2.txt
python3 -B verify_psi.py n12-bfs-center-s1.txt
```

生成命令（输出改为新路径，Clang模块缓存按项目约定设置）：

```sh
c++ -O3 -std=c++20 -pthread psi_search.cpp -o /path/to/temporary/psi_search
/path/to/temporary/psi_search 8 36 2 10000000 /path/to/new-n8.txt 1000 bfs
/path/to/temporary/psi_search 12 78 1 10000000 /path/to/new-n12.txt 1000 bfs
```

C++浮点搜索用于发现候选；精确证据是保存词的独立置换/峰判定及Fraction值。
同名resource文件是实际原生产资源回执，不是本轮轻量复核耗时。批次尚未全部完成。
