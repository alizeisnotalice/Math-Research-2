# R29 实际运行记录

本机 Python 3.15.0、macOS 27.2 arm64；纯整数 Walsh 与二元向量运算。最终源由独立审计 Agent 实际运行、正常退出，报告耗时约0.09360秒。主 Agent 阅读枚举闭包逻辑和结果，并运行归档验证器。无随机种子、无数值容差。

- 两块反例：k=2,3,4,5,6，各全部8相位，秩3、bent性及两块频率交换全部通过。
- 四块取等：k=3,4，各全部128相位，秩7，8种斜率置换固定像均通过。两类合计296相位。
- 四个有限集的全部置换及包含恒等置换的全部合法方向子空间：S={0,1,2}、{0,1,2,3}、{0,1,2,4}、{0,1,2,4,8}；置换数6、24、24、120；子空间数6、60、27、155；最大维数1、3、2、2。算法从零空间逐个加方向并检查整个新陪集合法，故不会漏掉该有限集合中的子空间。计数只指含恒等基点的族，不声称计数所有不同仿射陪集。

复现：

```sh
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round029/check_r29_counterexamples_and_permutation_families.py
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round029/verify_artifacts.py
```

脚本在自身目录写JSON。无本轮MacBook、GPU或MATLAB作业，未计入其他任务的占用。一般结论由解析归纳与循环证明承担；本轮未声称这些数学结论的外部首次性。继承来源（重新核对适用条件）及SHA256：

- `fourier-sparsity-dimension/proofs/002-affine-partitions-zh.md`: `57b256b3f50f7de127da6b9321038c67cfd3b88733edb9a9f2beed967d896a5e`
- `fourier-sparsity-dimension/proofs/006-disjoint-switch-closure-zh.md`: `2edb22747680e019e80d1551010d87fea82c9f6beb7aaf9ea74f754850abb984`
- `fourier-sparsity-dimension/proofs/007-minimum-switch-rank-zh.md`: `c897311431110b42cbbf79e537ef6d9783834d19a728f92aec8d46c55ddabb1a`
- `fourier-sparsity-dimension/proofs/017-double-block-rank-budget-zh.md`: `40717e76ad7bc9ee12dc35f7da536c53c20c9965bdf613ace4128aae18e2367f`
