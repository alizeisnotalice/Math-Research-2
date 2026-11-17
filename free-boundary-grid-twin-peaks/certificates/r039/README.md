# 第39轮证书

chain-obstruction：局部链引理、已否定的纯链覆盖候选及逐行反例。mixed-mechanisms：两项严格分解及row/comb的不同代价，混合全称界仍待证。warm-search：只增加保存词初态的搜索源及同态/坏输入门槛。finite-inputs：四个实际已完成链的续算初态、精确验证、资源和来源；不包含未完十亿步输出。

```sh
python3 -B free-boundary-grid-twin-peaks/certificates/r039/replay.py \
 --output /Volumes/data/CodexTemp/T2/a022-r039-your-replay
```

需要Python标准库、clang++、第36/37轮冻结源。输出目录全新且在data卷临时目录；编译缓存使用data卷共享clang缓存。重放9个根检查，门槛内部另记录18个进程；不启动长生产或GPU。
