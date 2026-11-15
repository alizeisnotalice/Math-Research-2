# 第8轮 GPU 探索：必要源码、原始记录与更正

状态：数值探路，不是原双峰条件律的概率估计，也不是严格完成数证书。严格 n8 结论来自独立 C++/Python 整数枚举。

## 算法与诊断

importance 两目录按当前外边界均匀选点生长，用边界基数之积作重要性权重；理论恒等式针对理想精确均匀随机数。GPU 实现使用有限精度对数权重及 xorshift32 取模，只有数值用途。11 个固定种子轨迹完成 CPU/GPU 对照。n4/n5 校准接近精确值；n12 权重集中、ESS 很低，不能外推。

smc 系列以边界基数作 potential、CPU systematic resampling、GPU 传播。跨批次必须以 Zhat 加权：sum(Zhat*p_hat)/sum(Zhat)，不能直接平均 p_hat。源含固定种子；不同种子不等于已证明独立或无重叠随机流。ESS 只是诊断，未认证混合或误差。继承 tag 的谱系严重退化时，粒子总数不是有效独立样本数。

`python3 -B summarize.py` 从原始记录重算 `summary-recomputed.json`。n7 原4点配置128批得到3.24518（精确3.2795）；n12同配置Z-ESS仅6.04/128、最大权重32.98%，该配置不可靠。边界延长L族 n7/n8 的32批估计3.3649/3.7707，相应精确值3.43617/4.05523；n10/n12 Z-ESS约4.4/5.3（共32批），不作阈值判定。有限批次偏离不等于已经证明系统偏差，也不能断言所有SMC算法失效。

## 必须保留的错误更正

lshape 第一批的实际外边界为7，历史输出 `scaled` 硬编码乘6。原始记录保留，p/logZ不因显示倍率而变；重算改乘7。其n7校准3.890208接近精确3.881127，不能沿用曾误称校准失败的判断。n8第一批修正值4.04571仍仅是数值候选。batch2已自动核验边界，重算4.61236；批次不稳定，不等同严格反例。这里的4点L与严格证书中的15点延长L不是同一输入。

硬件计时来自 Metal commandBuffer.gpuEndTime-gpuStartTime；wall字段是提交等待时间，不含全部CPU计算和初始化，因此不是完整程序墙钟。不能与第6轮旧 `gpu_seconds`（实际墙钟）直接合并。旧源码有些批次未显式检查commandBuffer.status；保存的CPU校准仅覆盖其明确测试，不能冒称全部GPU过程独立验证。

importance-extension/performance_during.txt 保存25个AGX运行采样，23个设备利用率至少94%，其中3个为100%。这是设备全局瞬时指标，不证明该任务独占或全程满载。

## 复现与资源

Swift 6.4.2 / Apple M6 / Metal 运行时编译。选择全新临时输出目录，避免历史文件被覆盖。对任一上述Swift源，可按以下模式编译：

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules SOURCE.swift -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/PROGRAM
```

各源的最后数行固定输入、粒子数和批次；SMC可传全新输出目录，importance版本按源码的命令行/当前目录输出约定运行。所有已归档批次已完成并退出；新一轮GPU批次另有PID，不计入本轮。源码、必要原始结果为自写实验产物，没有归档二进制或缓存。
