# 本机MPS深窗口96保存样本

m64/n104,d15,B8192,steps36000,restart6000,seed2026111704；六窗口各16候选，共96。worker实际结束，PID70700已消失，release记录清理调用及complete。父用Fraction重新计算每个保存的十进制样本、其所在窗口及前缀，96项均非精确反例；最小精确残差对应约0.00119764311239。数值优化阴性不代表有限排除；46,011前缀池的完备性继承已验证输入。

公共prefix-plan、基础worker、fractional_prefixes复用047-gpu-window15；Fraction验证器复用040-GPU-d14/verify_prefix_fraction.py。封装重算命令（在本目录）：

```
python audit_deep_windows.py --run . --worker worker.py --prefix-plan ../047-gpu-window15/prefix-plan.json --base-source ../047-gpu-window15/frozen-base-worker.py --verifier ../040-GPU-d14/verify_prefix_fraction.py --out recheck.json
```

恢复用state.pt约13MB为可再生临时检查点不入Git；最小复现源、输入绑定和所有96保存样本已封存。后续m80/n120/seed2026111705为独立运行批次，不计入本包。
