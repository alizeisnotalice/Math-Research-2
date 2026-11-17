# 本机d14 MPS检查点恢复

原任务随本机2026-11-17 05:01:02重启中断。只使用state.pt的step7000/epoch3（SHA c4f0f2f2558486bf2d9060253f1362ffe6066b7c40eb2a294784dbf772bd3652），舍弃晚于状态的step7100日志。新目录实际恢复运行；本目录不声称生产批次已经结束。

保存原执行字节不修改：resume_gpu_factor_search_prefix.py依赖的原工作器和前缀模块在frozen-worker中；portable_resume.py只把加载路径指向本目录副本。它保留原worker哈希校验。当前完整checkpoint保存在data上的任务目录，未将大体积可再生训练状态公开提交。

CPU预检核对源、全部参数、前缀、tensor形状/有限性/范围、Adam超参数与精确更新计数、完整重启链；六类语义负测通过。Adam参数和矩恢复自实际状态；旧RNG未保存，所以未来随机重启不能宣称逐比特延续。4小时到时原子保存partial状态并释放，只有12000步结束才可complete。

复现入口（准备旧检查点文件及其config/restart链后）：

```sh
python portable_resume.py --m 60 --n 100 --starts 32768 --steps 12000 \
 --fractional-index 14 --margin .001 --restart-every 2000 --seed 2026111615 \
 --resume-run <旧完整检查点目录> --out <新的data临时目录> --seconds 14400 \
 --validate-resume-only
```

去掉最后一个参数才运行GPU。有限浮点搜索不提供严格排除；原有实际源/状态预检回执、配置和恢复provenance已归档。完成批次仍需独立Fraction验证保存系数及确认进程退出。
