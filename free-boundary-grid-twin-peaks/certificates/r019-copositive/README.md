# 单调函数负相关的有限锥证书

证明 ../../proofs/monotone-root-copositive-finite.md 给出完整有限面枚举的充分性。所有计算使用Fraction；不依赖外部SDP收敛定理。

```sh
python3 exact_simplex.py ../r019-cut-full/n3-ie.json ../r019-cut-full/n4-cut-covariance.json ../r019-cut-full/n5-cut-covariance.json ../r019-cut-full/n6-cut-covariance.json ../r019-cut-d4/n7-d4-cut.json
```

共119个面驻点系统，实际都非奇异；n3..7得到严格负的最大值。四个独立小矩阵覆盖奇异平面、常数函数、边界/内部极值，结果见self-checks.json，输入同时归档。这个有限共正结论不是负半定结论，更不是对n一致的负间隙。

复现时请先确认 data 挂载，将 TMPDIR/TMP/TEMP 及编译输出设到 `/Volumes/data/CodexTemp/T2`；Swift/Clang 模块缓存为 `/Volumes/data/CodexTemp/cache`。环境：Python 3.15，Apple Clang 21，Swift 6.4.2/Metal。源码相对路径调整不改变数学核；原运行源定位保留在本轮报告，归档哈希见 manifest。有限核验不外推到所有 n。
