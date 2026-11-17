# MacBook 52次GPU中断前缀审计

24个种子目标；两完整目标各67,108,864 mask，第三目标连续前缀23,166,976/33,554,432。总收到157,384,704数值mask、19,212保存样本。所有未收到部分未知，原批次无终止/退出回执。2026-11-17恢复SSH后旧控制器和Wolfram内核均不在，但不能由进程消失补造成功回执。

主线程全量重新分配根组重构保存系数，41.836秒通过；GPU/CPU最大差7.0863e-5，根组乘积误差1.5749e-4。唯一近例为task0/mask525312的预设公平控制，四舍五入后的整数乘积已精确核对。六项语义负测试通过。所有浮点数据只作探索。

`stdout-*.jsonl.gz` 按完整行分片。复现时先把 inputs.json、search.wl、remote_gpu.py、payload.py 复制至 data 的独立临时目录，再执行：

```sh
python unpack_stdout.py --out <临时目录>/stdout.jsonl
python audit_gpu52_partial.py --root <临时目录> --out <临时目录>/audit.json
```

分片重建已由主线程实际执行，4,840行，完整SHA256为20ad4917f47d5aa744b22439fa14e32884d0422ad9edfe1ddb877772a3ccffa9。分析器记录所有原输入/源/流哈希，重新计算而不信任成功标记。保存输出和代码为本项目产生，不归档文献全文。
