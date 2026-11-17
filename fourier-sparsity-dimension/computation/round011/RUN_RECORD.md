# 第11轮实际核验记录

日期：2026-11-17。所有数学谓词均为整数/F₂运算，计时使用浮点。

## 主Agent实际执行

- `verify_rank_and_code_bounds.py`：128输入选择树的完整FWHT，谱16点、线性秩7、仿射秩6；K₄(6)=K₄(10)=−20及第二字符矩；真实六维原子码全部64码字、4096对自正交检查及16个重量六字的Gram矩阵；模选择器8个合法相位与一个非法乘积负对照。0.030923秒。
- `find_dim7_local_model.py`：枚举16原子七维码的128字，并完整遍历65536个候选相位，逐一测试至多四个bent条件。无兼容Boolean行。此模型不是实际残余；它用于区分码重预算与Boolean合法性。
- `verify_graph_spectra.py`：每份证书仅作一次长度32768的整数FWHT直接验证输出函数，另作单系数变号负对照；20份全部通过。完整枚举支撑平移，10份周期空间为零、10份一维。0.849196秒。输出含每份完整支撑与系数分子，任一零周期记录就是可独立复验的有限反例。
- `verify_receipts.py`：绑定实际源文件、历史最小输入、修正后的本机/远程输出以及独立输出域算法，核验全部20条跨机结果相等，并与当前仓库中的源证书哈希比对。

## 两机逐切片核验

`audit_common_affine_derivative_periods.py`读取`period-input.json`，恢复七个全局仿射基及八个相对生成元。完整枚举64个方向，检查基函数导数常值及图平移判据；对10个零周期例子，各重算32768个六变量整数Walsh变换。

修正后的同一源程序SHA-256：5b4490874a1b64be61679d833fec8a07a89152cf191f1cce3b16ba38598c6699。本机Python3.15.0实际耗时54.977288秒；MacBook Python3.9.6实际耗时73.998046秒，主机名zzhMateBook-Pro-S.local。远程由计算子Agent用一次SSH stdin bootstrap及nice 10运行；其执行回执为退出0。主Agent实际取回完整JSON，核对源/输入哈希及20条结果，并检查远程任务路径无残留进程，见`remote-fetch-and-exit.json`。

MacBook本轮为单进程独立复验，不称CPU满载；第10轮的15/30路负载属于上一轮历史。未启动MATLAB、Mathematica或GPU上下文。全部本轮计算已结束。

初版周期脚本曾把64位函数掩码误用六位秩检查，已在输出前修正；另一版未把失败字段传播为非零退出。归档源程序使用64位掩码秩、失败assert及Python3.9兼容位计数，归档两机结果均来自该最终源。旧62.432秒输出及中间版本没有作为当前验证归档。

## 输入与复现

`period-input.json`是第6轮历史输入的原字节副本，保留其历史计时与来源字段；本轮只读取seed permutation、h_bits、relative_generators和来源定位。不得将文件内旧计时视为第11轮计算。当前结果按完整源哈希记录，所有20份相位与基可由此输入重建。

在仓库根目录执行；输出请改到data盘临时目录，避免覆盖已归档回执：

```sh
python3 fourier-sparsity-dimension/computation/round011/verify_rank_and_code_bounds.py --out /Volumes/data/CodexTemp/T2/a28/r011/rank-replay.json
python3 fourier-sparsity-dimension/computation/round011/find_dim7_local_model.py --out /Volumes/data/CodexTemp/T2/a28/r011/model-replay.json
python3 fourier-sparsity-dimension/computation/round011/audit_common_affine_derivative_periods.py --audit fourier-sparsity-dimension/computation/round011/period-input.json --out /Volumes/data/CodexTemp/T2/a28/r011/period-replay.json
python3 fourier-sparsity-dimension/computation/round011/verify_graph_spectra.py --input fourier-sparsity-dimension/computation/round011/period-input.json --out /Volumes/data/CodexTemp/T2/a28/r011/graph-replay.json
python3 fourier-sparsity-dimension/computation/round011/verify_receipts.py --out /Volumes/data/CodexTemp/T2/a28/r011/receipt-replay.json
```

周期/图程序兼容Python3.9；加权码模型程序使用int.bit_count，要求Python≥3.10，本次实际版本为3.15.0。一般环境的维数定理与迭代定理由证明承担，不由这些有限运行承担。文件哈希见manifest.json。
