# r11 GPU：高状态随机流、逐层矩与组内归约核验

Apple M6 / Metal、Swift6.4.2。生产使用冻结文件 source/gpu_tau_xoshiro256ss_production_r11.swift；source/gpu_tau_xoshiro256ss.swift 添加独立 audit 入口，Metal 内核相同。后者 n7/n8 每τ检查首个256粒子组：存活数精确，存活组的六个log矩容差2e−4；全死组单独核对存活数，零贡献由完整批次原始矩另核验。每τ单轨迹的原始64位随机字、顶点顺序、峰距均与CPU精确一致，log权重容差1e−4。

算法使用 xoshiro256** 1.0 + SplitMix64 填充状态，拒绝法选集合索引。出处：[作者参考源码](https://prng.di.unimi.it/xoshiro256starstar.c)。主 Agent 核对 next 转移、旋转常数7/45及作者的SplitMix播种建议；多个seed不等于数学独立性或无重叠证明。n7两生产域各33554432轨迹；n8八生产域各50331648轨迹。原始矩使用未除n²的e2，报告Q与其经验方差须相应除n²与n⁴。

生产与audit分开，失败debug/smoke不计入。元数据、batch、stratum及summary全部保留；零质量层null表示log0，样本分母仍包含死路径。实际结果由课题 verify_r011_gpu.py 独立汇总，不能作为原问题的概率认证。n8 pooled Q≈.334175405，简单重复均值≈.334192687，两者不同。

每rep的gpu_seconds是硬件执行区间和，wall_seconds是commit+wait之和，不含编译、CPU核验/汇总和批间空档；没有单独记录完整进程墙钟。n8每rep硬件14.3–15.8秒，提交等待14.85–17.39秒。该批未连续记录AGX；主Agent在第4个rep运行时观察整机96%单点，不能归因于单一任务。所有这些生产与audit进程已退出。

复现时在此目录运行（输出目录必须尚不存在，编译件放CodexTemp）：

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules source/gpu_tau_xoshiro256ss_production_r11.swift -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r11-gpu-sis
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r11-gpu-sis /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/reproduce-r11-n8 metadata/n8-adaptive-allocation.json a022f80000000000
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules source/gpu_tau_xoshiro256ss.swift -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r11-gpu-audit
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r11-gpu-audit /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/reproduce-r11-audit metadata/n8-adaptive-allocation.json a022f80000000000 audit
```

其他生产seed取各run首行metadata，固定allocation保存于metadata。原结果中的绝对路径只记录历史输入，复现可用本包相对metadata路径。二进制、缓存和失败调试输出未提交。
