# 第12轮实际核验记录

日期：2026-11-17。主机为本机Mac mini。归档输出均绑定实际源SHA-256；没有远程成功运行回执。

## 已执行的最终程序

- `verify_module_partitions.py`：Python3.15.0，F₂真值掩码与乘法余式核，完整检查Q=2,4,8的2+15+4140=4157个分割，0.152796秒。直接计算最大保留模块维数，再与同色内部节点数及预算比较；没有把节点计数本身当作独立核验。
- `recover_weight_lp_vertex.py`：Python3.13.18、SciPy1.18.1、SymPy1.14.0，种子2026111712，首个候选成功，0.034288秒。HiGHS浮点候选后，取13个独立活跃等式作精确有理解并直接检查枚举约束。活跃检测阈值1e-5、零坐标阈值1e-7只控制候选选择；最终判定用Fraction，零容差。
- `verify_residual_and_lp.py`：Python3.15.0，0.057619秒。独立组装全部65个MacWilliams重量变换及正负包含条件；同时核验31个包含1_R的五维真实残余子码，原子数为12的30个、16的1个。完整参数行均Boolean，原子变换支撑、Plancherel与重量条件通过。精确确认61/768−79/1024=7/3072及四个尺度的K₄公式。

三个最终源程序和对应JSON均在本目录；数学判据用整数/F₂/Fraction，计时用浮点。模块验证覆盖一个基代数原子的全部小分割，任意尺度、多原子的结论由证明承担。LP计数非整数，仅是连续松弛可行证书。

## 探索纠错与退出

初次浮点D扫描因尺度问题不作结论。一个早期活跃集恢复误把B≥|K|写成B≥−|K|，主Agent独立全重量检查发现负B₄，已拒绝该证书。最终程序修正为正确绝对值方向，且独立检查所有65个重量。旧候选与中间输出未作为数学证据归档。

另一个完整精确单纯形探索运行8分29秒仍未结束，主Agent核对本任务命令后发送TERM，ps退出码1确认进程已消失；这不意味着LP不可行。回执见`exact-simplex-termination.json`。所有本轮计算已结束；没有启动MATLAB、Mathematica或GPU上下文。

## MacBook

已有SSH别名与重新发现的IPv4、IPv6地址均连接超时，尚未进入认证。由远程Agent再次有界核查，结果相同。没有启动新远程计算，没有读取到CPU负载，未改SSH或系统设置。摘要见`remote-unavailable.json`；详细网络诊断保留在本机临时目录，不提交设备地址配置。第6至10轮的高负载记录属于历史，不能当成本轮满载证据。

## 复现

从仓库根目录执行；输出使用data盘临时目录，避免覆盖归档回执。首个命令需现有SciPy/SymPy环境，其余仅标准库。

```sh
/Users/zzh/.venvs/main/bin/python fourier-sparsity-dimension/computation/round012/recover_weight_lp_vertex.py --out /Volumes/data/CodexTemp/T2/a28/r012/lp-replay.json
python3 fourier-sparsity-dimension/computation/round012/verify_module_partitions.py --out /Volumes/data/CodexTemp/T2/a28/r012/module-replay.json
python3 fourier-sparsity-dimension/computation/round012/verify_residual_and_lp.py --certificate /Volumes/data/CodexTemp/T2/a28/r012/lp-replay.json --out /Volumes/data/CodexTemp/T2/a28/r012/residual-replay.json
```

输出哈希绑定与交付文件哈希见`manifest.json`。外部论文只记录来源、定位和哈希，未提交全文。
