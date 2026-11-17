# R28 实际运行记录

独立审计 Agent 在本机 Python 3.15.0、macOS 27.2 arm64 实际运行最终源，退出正常；最终报告计时约0.01998秒。整数 GF(2) 运算，无浮点容差。主 Agent 阅读最终源及输出，并运行归档核验器。

- K₄：ρ=4，N=r=10，W维6、q=64；立方体图：ρ=8，N=r=20，W维12、q=4096。逐行重量4、秩、实际点标签互异及标签异或和0均通过。采用稀疏标签，不枚举指数大小的所有相位。
- 80个W=F₂⁵四点生成模型，随机种子20261117、至多10个输入生成元。先选独立的原重量4生成元，再完整枚举每个小对偶码。最高秩10，最高N=27、r=17；全部列界通过。这不是全部32点关系空间枚举。
- 归档前修正了初稿把依赖生成元展开长度记作码大小的问题；最终源按独立基展开并核验D₄秩与D秩相同。只有修正后实际输出作为本轮证据。

复现：

```sh
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round028/check_d4_residual_bound.py
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round028/verify_artifacts.py
```

程序在自身目录写结果 JSON。无GPU、MacBook或MATLAB任务；不把别的课题的CPU负载归入本轮。一般结论由解析证明承担，源哈希只是文件绑定。继承定理此次重新阅读了适用假设，定位与SHA256如下：

- `fourier-sparsity-dimension/proofs/007-minimum-switch-rank-zh.md`: `c897311431110b42cbbf79e537ef6d9783834d19a728f92aec8d46c55ddabb1a`
- `fourier-sparsity-dimension/proofs/017-double-block-rank-budget-zh.md`: `40717e76ad7bc9ee12dc35f7da536c53c20c9965bdf613ace4128aae18e2367f`
- `fourier-sparsity-dimension/proofs/026-four-relation-closure-zh.md`: `eb1139bf21ee3e724b06da34e95aedfc5a36eb15db9905f41e57f839404d626f`
- `fourier-sparsity-dimension/proofs/027-four-relation-pair-budget-zh.md`: `976f4ba6a3c92f28da908d948dc12a0ba7b24433450b2cfaed6e7c6a5c109266`
