# r020：距离—低峰出生秩联合账本

`n3_direct.py` 全枚举 9! 个标号；`distance_birth_dp.cpp` 用确定性子集 DP 计算
\((d,\tau)\) 联合整数表。`verify_and_summarize.py` 交叉两法，并对照既有精确
\(A_2\)、峰距直方图，以及 n=3 已归档的出生秩直方图；它输出每个距离层的
\(\mathbb E\tau\)、\(\mathbb P(\tau=2)\)、\(G=N-\tau\) 质量和精确尾数。

本轮只执行 n=3、n=4，未启动 n=5。没有随机性或种子；使用一个 CPU 进程、
小于 16 MiB 的 n=4 DP 工作数组。

```sh
TASK_TMP=/Volumes/data/CodexTemp/T2/a022-distance-birth-replay
mkdir -p "$TASK_TMP"
c++ -O3 -std=c++20 -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang-modules distance_birth_dp.cpp -o "$TASK_TMP/distance_birth_dp"
python3 n3_direct.py > "$TASK_TMP/n3-direct.json"
"$TASK_TMP/distance_birth_dp" 3 > "$TASK_TMP/n3-dp.json"
"$TASK_TMP/distance_birth_dp" 4 > "$TASK_TMP/n4-dp.json"
python3 verify_and_summarize.py
```

结果只说明有限 n=3,4；特别是短距层中 \(G\) 的有限分布不能外推为一般尺度结论。
