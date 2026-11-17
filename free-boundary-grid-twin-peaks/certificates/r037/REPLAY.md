在仓库任意目录执行（输出必须为data卷下全新目录）：

```sh
python3 -B free-boundary-grid-twin-peaks/certificates/r037/replay.py \
 --output /Volumes/data/CodexTemp/T2/a022-r037-your-replay
```

需要本机clang++、Python标准库、bash和第36轮冻结搜索源码。缓存显式使用data卷共享clang缓存。此命令只做有限核验，不启动生产扫描或GPU任务。低支审计还会在其证书目录写一份可再生JSON；其完整内容同时保存在stdout。
