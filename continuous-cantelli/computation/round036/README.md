# 第36轮复现

quartic_jet_check.py 用 SymPy 和 mpmath 核对 Hermite 因子及误差的幂次；quartic_lagrange_check.py 独立检验纯四次模型。用 Python 运行各脚本即可。

从各台机器对应的子目录临时副本运行。两台机器的版本与第35轮相同，计算参数由下列命令和程序内确定性网格指定。

* metal-local 的 Apple 原生程序使用当前目录中的 quartic.metal。将 quartic.swift 以 swiftc -O 及已有共享模块缓存编译为 quartic，再执行 `./quartic 90 96 main 0`。verify.py 的两个参数依次为 GPU 输入 JSON 及 CPU 输出 JSON。exact_witness.py 复核在同一测试中跨过阈值的有理证书。
* cpu-local 在共同窗口内复用第35轮 high_jets.py 的12个进程。synchronized-local-summary.json 列出 PID 和重叠采样。heat_potential.py 另做了慢平坦轮廓的小规模对照，参数和精度见 smoke.json。
* metal-remote 以 swiftc -O 编译 scan.swift，当前目录须有 dual.metal。正式运行的参数网格见 sync.json。用 `python check_cpu.py sync.json sync.cpu.json` 复核选定参数。rational_gaussian_well.py 按有理分数方式验证固定尾井。
* cpu-remote 设 BLAS 为单线程，并行 15 份 `python finite_packet.py i --exponent 16 --order 32 --seconds 89 --max-repeats 160 --out shard-i.json`，以 0 至 14 替换 i。exact_n4.py 执行精确 Sturm 根计数，high_precision.py 独立复查选定积分。两台设备时钟差为事后往返测得，见 clock-comparison.json。
