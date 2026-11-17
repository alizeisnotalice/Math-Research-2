# MacBook Native Metal 数值基准

固定源码 NativeScore.swift、输入 task2-input.json 对应实际 MacBook Apple M5 Pro、macOS27.3、Swift6.4.2运行。remote-six.json保存六个明确mask，remote-1m.json保存重复范围[0,1048576)的GPU top4。remote-1m.time为/usr/bin/time实测0.49秒，约214万mask/秒；只是一次短基准，不计任何新增覆盖，也不作为长期吞吐承诺。

主线程实际运行 parent_check_native.py，从root-task2.json中的原始数值根组直接乘积重建九个非控制保存样本（六mask中五个加大批top4），对比GPU分数，最大误差0.000512638，小于0.001回归门槛；卷积残差另逐项记录。这个检查不证明整个大批top4正确、不证明精确代数分解。mask0为公平控制，单列null不进入排名。搜索器与候选验证器分离，原来的verify_cpu.py也保留用于追踪实际初审；它引用原本机输入路径，父新验证器可从本包直接复现。

remote-observation-20261117T130911.txt实际采样Device Utilization99%，同时v3 Wolfram任务存在，因此只归为共享设备占用，不能归为native独占或持续满载。remote-postbench-exit-concise.txt记录native_score PID已不存在；正式v3任务继续。

示例复现（复制本包到有写权限的临时目录，在Mac上用已有Swift编译；不要在活动计算旁盲目重复大批）：

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift NativeScore.swift -o native_score
./native_score --input task2-input.json --masks 0,1,17,1234,23166975,23166976 --out six-new.json
./native_score --input task2-input.json --range 0,1048576 --out million-new.json
/Users/zzh/.venvs/main/bin/python parent_check_native.py
```

远端使用其已存在的CodexRemoteMath缓存目录，不能假设MacBook挂载/Volumes/data。本包原型对nfreq≤64、ncoef≤53有限制；不是任意输入的最终生产调度器。下一正式程序仍在开发，尚未接管v3。公开归档只选本题必要观察，未收录完整远端进程列表。
