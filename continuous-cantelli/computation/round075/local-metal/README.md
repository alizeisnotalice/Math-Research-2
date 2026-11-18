# R75 本机 Metal 实测

沿用第68轮的 Cantor 最小集圆周高斯模态着色器。本轮为
64 点 Gauss–Legendre，计数区间
`[536,870,912, 920,813,568)` 与前次的48点计数段不交。

- 本机 Apple M6 实际 UTC 为 2026-11-18T04:01:00.997Z 至 2026-11-18T04:03:01.002Z。
- 墙钟 120.005728 秒，Metal 命令区间 101.018386 秒，提交和等待 118.850985 秒。
- 共 383,942,656 组模态和参数。PID 28917 退出码 0，已确认退出，无非有限结果。
- 222 次整机 GPU 采样均值 100.000%，范围 100–100%。启动前为 99%，结束后为 100%。整机读数包含其他作业。
- 10 个分层代表经 CPU 独立积分复查。总模态、粗层、细层的最大绝对误差分别为 8.21e-07、8.5e-07、3.46e-08。这项核验覆盖所列样本。

## 复算

```sh
CANTOR_START_INDEX=536870912 ./run.sh 120 64 main
PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python verify.py cantor-main.json cpu-main.json --precise
/Users/zzh/.venvs/main/bin/python summarize.py
```
