# n≤11 CPU原律根矩

__uint128_t邻接/占据，121顶点前缀；线程局部拒绝取样、整数CDF。每个非死亡终端样本额外用直接邻点rank比较器检查恰两峰，并断言整数F−2G=4d²。F为两根双倍中心坐标平方和，G为双倍坐标点积；V/n²=ΣF/(8n²ΣZ)，C/n²=ΣG/(4n²ΣZ)。

```sh
clang++ -O3 -std=c++17 -pthread covsmc.cpp -o /Volumes/data/CodexTemp/T2/a022-cov
/Volumes/data/CodexTemp/T2/a022-cov 10 4096 32 16 8 /Volumes/data/CodexTemp/T2/a022-cov-n10.csv 0xA0222400
python3 summarize_covariance.py n10-p4096-r32-l16.csv
```

regression.json记录n3/n8全部秩、L0/4/16与前轮核同seed的所有旧字段逐字符串一致。n10小粒子ASan/UBSan测试通过。完整n10批3168记录，272.00秒、RSS19,300,352 B，退出后归档。字段采用双精度计数权重，后聚合Decimal80不会把双精度生成数据变为精确数。协方差以D4人口均值（几何中心）为准，不再减去随机样本均值。随后完整n8(P16384,R32,4线程)及n11(P8192,R32,8线程)批也已完成，分别2016和3840记录；三批都按原CSV直接Decimal80复算，资源与metadata一并归档。进一步n11 P16384,R64的独立域批次另放r020临时目录，尚未完成，不在本证书内。

复现时请先确认 data 挂载，将 TMPDIR/TMP/TEMP 及编译输出设到 `/Volumes/data/CodexTemp/T2`；Swift/Clang 模块缓存为 `/Volumes/data/CodexTemp/cache`。环境：Python 3.15，Apple Clang 21，Swift 6.4.2/Metal。源码相对路径调整不改变数学核；原运行源定位保留在本轮报告，归档哈希见 manifest。有限核验不外推到所有 n。
