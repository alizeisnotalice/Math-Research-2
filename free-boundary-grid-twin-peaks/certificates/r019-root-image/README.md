# 可拆根像 n5 与带纤维预算的散布支

完整推导在 ../../proofs/root-image-mean-and-scattered-insertion.md。以下重算每个首次合并状态的完成数与根邻接权重。

```sh
clang++ -O3 -std=c++20 -pthread root_image_n5.cpp -o /Volumes/data/CodexTemp/T2/a022-image
/Volumes/data/CodexTemp/T2/a022-image 5 /Volumes/data/CodexTemp/T2/a022-image-n5.json 4
python3 verify_r019.py
python3 audit_scattered_subsets_n3.py
python3 audit_scattered_lobes_n3.py
```

n5的I/A1≈28.4273只是有限精确比。I≥cNA1至多重得已有线性R下界，不能替代exp(c n log n)级增益。散布支映射严格满足每个输出总权≤1，但n3的平均增益小于原连续块映射；没有证伪原题。

复现时请先确认 data 挂载，将 TMPDIR/TMP/TEMP 及编译输出设到 `/Volumes/data/CodexTemp/T2`；Swift/Clang 模块缓存为 `/Volumes/data/CodexTemp/cache`。环境：Python 3.15，Apple Clang 21，Swift 6.4.2/Metal。源码相对路径调整不改变数学核；原运行源定位保留在本轮报告，归档哈希见 manifest。有限核验不外推到所有 n。
