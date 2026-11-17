# R31 实际运行记录

最终源由主 Agent 在归档目录实际运行：Python 3.15.0、macOS 27.2 arm64，exit 0，JSON计时约1.44633秒。纯整数FWHT及GF(2)；无随机种子、无浮点容差。

八变量模型基相位及四个秩4图块按第18轮公式构造，并提升到新增v=0层。全局旧仿射方向与块指标选出9维基，完整检查512个相位。

- 四块各32点、互不相交，全部方向支撑于128点R；q=16、M=8。
- 全部相位bent；每块的两个γ片各16点、非零部分Walsh值±16；八片互斥，始终分割128点对偶u频率为0的超平面。
- 各片的方向保持固定；按其法向记录标号平移，逐个参数直接核对其等于生成元记录的异或。联合平移秩5，核恰为四块常数空间。
- 对每块各方向的二次系数逐项核对：仅第七个旧全局x₅方向非零，其系数向量恰为Q=x₁x₂+x₃x₄，其余方向为零。模⟨pol Q⟩冻结秩为0，τ为该参数位。初稿只查一维变化，归档前加强为这个精确向量检查，并实际重跑最终源。

复现：

```sh
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round031/check_r31_scalar_rank4_absorption.py
/Users/zzh/.local/bin/python3 fourier-sparsity-dimension/computation/round031/verify_artifacts.py
```

脚本在自身目录写JSON。本轮没有MacBook、GPU或MATLAB作业，计算进程正常退出；未复述历史远程运行作为本轮证据。有限实例验证h=1及非平凡τ，不覆盖h=2,3的所有图族，不证明新容量界取等。一般结论由解析证明承担。继承来源与SHA256：

- `fourier-sparsity-dimension/proofs/018-quadratic-graph-blocks-zh.md`: `c508732d4c7fa5655c0bc4afbe7387e4842afe0bdc7ab7502bd2586f0b953b49`
- `fourier-sparsity-dimension/proofs/019-projectable-double-blocks-zh.md`: `b66a34ddd5cda30c8cb7aa2ca57030b6b44b0f80a202aee773d34e8594650495`
- `fourier-sparsity-dimension/proofs/030-moving-coset-local-budget-zh.md`: `bdd6efbb645683a4eeea301ac77de7425f4be1259eb2e376c8af5b164ea73fb2`
