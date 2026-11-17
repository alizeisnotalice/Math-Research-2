# 复现

需要 Python3 标准库与支持 C++20 的 Clang。所有可再生输出写入指定的新目录；
`--module-cache` 指向共享 Clang 模块缓存。示例从本证书目录运行：

```sh
python3 -B replay.py --out /Volumes/data/CodexTemp/T2/a022-r035-replay-new \
  --module-cache /Volumes/data/CodexTemp/cache/clang
```

该命令实际重放四项有理/枚举审计、保存的 n7 CRT 和残差，重新编译空间核，运行
n3全部9孔与n5六孔类共60个小模任务，并与独立完整联合表比较全部25字段。
不执行 n7 生产，也不运行 n11 或任何 GPU 批次。每个子进程 wait/reap 后检查退出。
保存的原生产耗时不算本次复现耗时；临时工作编号 r037–r041 不代表公开轮号。
