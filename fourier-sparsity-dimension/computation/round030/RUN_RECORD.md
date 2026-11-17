# R30 实际运行记录

本机 Python 3.15.0、macOS 27.2 arm64，纯整数FWHT/GF(2)。主 Agent 实际运行归档目录中的最终源，exit 0，结果计时约0.07825秒。无随机种子、无浮点容差。

- 四块例：全部128个相位bent；空间维7，支持32点，四个8点仿射块确有两种不同方向。全部相位的四部分谱均0或±8、互斥，合成固定32点D_R（对偶u坐标为0）。平移方向秩3，核恰是4维块常数空间。
- 三块例：全部16个相位bent；空间维4，支持24点，三块方向不全相同。总部分谱固定24点，因基数而非仿射；平移秩1，核恰是3维块常数空间。
- 两例还直接检查块指标的方向基坐标、块两两不交、每个方向支撑在R内及块的实际仿射重构。

复现：

```sh
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round030/check_r30_small_examples.py
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round030/verify_artifacts.py
```

源自独立审计 Agent 初稿，主 Agent 归档前发现并纠正了参数斜率未乘参数位、误把32点部分谱当作64点全谱、用物理块不同代替方向不同的实现问题；整合修正时又纠正了一处已返回方向空间却重复平移的判断。只有上列最终源实际通过的输出作为计算证据，初稿失败不计入通过范围。

数学主稿由独立理论 Agent 审查；混合q/2q块的一般界由解析冻结与核分账证明承担，本轮有限例子仅覆盖最小块。无MacBook、GPU或MATLAB任务，进程已正常退出。另有理论 Agent 的小维边界扫描未纳入本轮核心证据。继承来源及哈希（已复核适用条件）：

- `fourier-sparsity-dimension/proofs/002-affine-partitions-zh.md`: `57b256b3f50f7de127da6b9321038c67cfd3b88733edb9a9f2beed967d896a5e`
- `fourier-sparsity-dimension/proofs/017-double-block-rank-budget-zh.md`: `40717e76ad7bc9ee12dc35f7da536c53c20c9965bdf613ace4128aae18e2367f`
- `fourier-sparsity-dimension/proofs/029-finite-permutation-residual-zh.md`: `0b7770266185fd192bdde8172446e554d4be74aeb4c752d679733fdc260b1625`
