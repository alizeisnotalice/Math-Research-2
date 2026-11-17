# R27 实际运行记录

本机 macOS 27.2 arm64，Python 3.15.0，纯 Python 整数算术，无容差、无随机种子。两份脚本由独立审计 Agent 实际运行，结果随最终源归档；源内容未改，主 Agent 逐段检查源及输出。未把源哈希核验当作数学证明。

1. `check_fourier_sparsity_le4.py`：全部65536真值表，1152个稀疏≤4者结构检查无失败，计时约0.227秒。
2. `check_twin_closure_and_pairings.py`：k=3,4,5 的pair码全部8、128、32768个对偶字；两个k=4小模型；RM1(3)对照。另检查896个bent种子×9种配对=8064例，每例4块、64频率，部分谱幅度及逐频分割全部通过，计时约7.801秒。这是有限模型证据，不是所有种子/共同家族枚举。

复现（脚本在自身目录写 JSON）：

```sh
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round027/check_fourier_sparsity_le4.py
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round027/check_twin_closure_and_pairings.py
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round027/verify_artifacts.py
```

两项计算均正常退出；本轮没有 MacBook、GPU 或 MATLAB 作业，没有持续 A28 CPU 满载声称。完整一般预算依赖本项目第26轮闭包及第17轮部分仿射块定理；本轮重新阅读相关假设和账本，并独立审查商空间残余传递。源定位与 SHA256：

- `fourier-sparsity-dimension/proofs/026-four-relation-closure-zh.md`: `eb1139bf21ee3e724b06da34e95aedfc5a36eb15db9905f41e57f839404d626f`
- `fourier-sparsity-dimension/proofs/017-double-block-rank-budget-zh.md`: `40717e76ad7bc9ee12dc35f7da536c53c20c9965bdf613ace4128aae18e2367f`
