# A28 第10轮实际运行记录

日期：2026-11-17。全部数学谓词使用整数或F₂运算，浮点只用于计时及CPU占用率。下述历史记录指本轮实际执行，不是待执行计划。

## 程序和核验

- `verify_seed_independent.py`：独立Python有限域乘法、迹与整数FWHT；十四变量种子真值表逐位与C++一致，耗时0.709091秒。模多项式0x4021，α=0x13f2；Walsh直方图−128:8128、128:8256，ANF次数4。**没有执行非弱正规性穷举**。
- `cddl_fact22_lift_verify.cpp`：所有前缀布尔真值表的三角位移，逐w任意仿射x相位。每个样本计算完整真值表和整数FWHT，逐系数与独立闭式公式比较，并检查单点翻号负对照。随机种子、全部样本参数或其确定性生成输入随回执保存。
- 第一批MacBook14进程各运行ℓ=0…4、每层64样本：4480样本，直接检查5,005,901,824个系数；全部退出0。峰值941.7% CPU，两个整机采样分别12.18%、31.46% idle，不能称为满载。本机同一源程序执行其中320样本，逐项参数及结果与第0远程作业一致，耗时5.116854秒。
- `verify_residual_atom_constructions.py`：实际核验t4残余嵌入的三个补维尺度、五变量非仿射谱、六变量真实部分覆盖的全部1024个切片；0.216379秒。
- `verify_dim6_residual.py`：实际核验13维空间全部8192切片、64个残余部分Walsh、完整1395个线性三空间及11160个仿射三平坦。得到40个各自可对整个V作一维扩张的最小翻号；0.620860秒。完整基掩码和40个见证随JSON保存。
- `classify_fixed_atom_quadratic_extensions.py`：主Agent整合后实际重跑32商类×64平移，枚举全部可含子空间。六个好类、五条非零线、无二维扩张；每个非零好陪集均含四个重量8字。只覆盖固定原子二次模型。

## 扩大计算与CPU实测

用户要求充分调用MacBook，故另以新随机种子执行扩大核验和调度对照。扩大的原始逐样本JSON保留在远程任务专属目录；公开回执保存源程序、编译器、运行命令、确定性种子、原始输出哈希、全部检查汇总、PID、连续CPU采样及进程退出码。最小复现输入是源程序及CLI参数，不依赖缓存文件。

第一扩大批次15进程各每层1024样本，普通调度。执行中只对本任务PID调用`taskpolicy -B -p PID`，结果保存在`remote-scheduling-adjustment.json`。第二批15进程各每层512样本，第三批30进程各每层256样本；后两批采用`pthread_set_qos_class_self_np(QOS_CLASS_USER_INITIATED,0)`，除调度调用之外数学源代码相同。

精确汇总见`expanded-receipt-verification.json`。它核对源程序与调度脚本哈希、任务数、种子、返回码、逐作业成功谓词和系数数量；有限扩大抽样不替代正文对全部参数的证明。CPU满载仅以实际采样判断，不把工作进程数当成利用率。

MacBook型号及环境由回执实测：15逻辑CPU，Python3.9.6；本机Python3.15.0。编译器完整版本与命令见相应execution JSON。本轮没有新启动MATLAB、Mathematica或GPU上下文；这里的精确整数FWHT使用CPU。

## 复现

Python核验脚本要求Python≥3.10（本次实际使用3.15.0）。在仓库根目录运行，输出和编译文件放data盘的新目录；下例目录须预先存在。Python脚本输出参数可自行换到该目录，避免覆盖已归档回执。

```sh
clang++ -O3 -std=c++20 -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang fourier-sparsity-dimension/computation/round010/cddl_fact22_lift_verify.cpp -o /Volumes/data/CodexTemp/T2/a28/r010/replay
/Volumes/data/CodexTemp/T2/a28/r010/replay --max-l 4 --samples 64 --seed 202611171000 --out /Volumes/data/CodexTemp/T2/a28/r010/replay.json
python3 fourier-sparsity-dimension/computation/round010/verify_dim6_residual.py --out /Volumes/data/CodexTemp/T2/a28/r010/dim6.json
python3 fourier-sparsity-dimension/computation/round010/classify_fixed_atom_quadratic_extensions.py --out /Volumes/data/CodexTemp/T2/a28/r010/extensions.json
python3 fourier-sparsity-dimension/computation/round010/verify_lift_receipts.py --out /Volumes/data/CodexTemp/T2/a28/r010/receipt-check.json
python3 fourier-sparsity-dimension/computation/round010/verify_expanded_receipts.py --out /Volumes/data/CodexTemp/T2/a28/r010/expanded-check.json
```

旧探索程序中曾遗漏三角函数的非线性前缀项、且总成功标志仅包含种子检查；本轮归档C++已修正，所有归档执行均使用修正版。D6早期探索中的全局相位与零支撑相位混淆也已修正；归档程序直接核验真实零支撑掩码。被撤回的探索输出未当作证据。

来源文献的非弱正规种子性质保留外部计算认证依赖。计算阴性、线程调度和CPU实测均不承担一般A28结论。最终文件哈希见manifest.json；外部论文全文及可再生编译缓存不提交。


## 扩大批次最终实测

|批次|进程数|样本数|直接Walsh系数数|总墙钟秒|任务CPU峰值|整机最小idle|
|---|---:|---:|---:|---:|---:|---:|
|remote-full-cpu-execution.json|15|76800|85815459840|199.534|954.9%|9.26%|
|remote-qos-cpu-execution.json|15|38400|42907729920|86.779|1245.7%|9.24%|
|remote-qos-30-execution.json|30|38400|42907729920|87.926|1145.4%|16.73%|

三个扩大批次合计153600样本、171630919680个直接Walsh系数；连同首批总计158080样本、176636821504个系数。三个扩大批次均没有0% idle采样，不能称持续满载。15路USER_INITIATED的整机最高采样占用90.76%；30路并发没有进一步提升。本轮观察到其他计算进程同机运行，不据此直接断定空闲率的原因。只调整本任务进程/线程的调度，没有修改系统电源配置。

全部作业与SSH退出码0；最终对四个本轮专属执行路径的pgrep返回1且无输出，见process-exit-check.json。本次启动的计算已经结束，没有终止其他任务。
