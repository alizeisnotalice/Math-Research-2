# 第48轮复现说明

从本目录运行。使用 data 卷独立空目录保存二进制和重放结果。实际环境是 Python 3.15.0、Apple clang 21.0.0。整数账本与采样状态分开存放。

## 首次合并间距和路径

将以下命令中的 SCRATCH 指向数据卷空目录。

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Werror merge-repair/gap_dp.cpp \
  -o "$SCRATCH/gap_dp"
python3 -B merge-repair/audit.py "$SCRATCH/gap_dp" \
  | cmp - merge-repair/verification.json
python3 -B early-fibres/audit_early_fibres.py \
  | cmp - early-fibres/audit_early_fibres.stdout
python3 -B early-fibres/audit_thread_entropy.py \
  | cmp - early-fibres/audit_thread_entropy.stdout
```

间距 DP 按合并前的无边双分量、合并点、根的有序指派、低支的外接触目标集分解。高支、低支各用有根连通增长数，合并前交错数由最高峰固定词首确定。后缀 DP 记录第一次命中外接触集合的时间；命中后改用全部连通增长的计数。三阶原词穷举与上述程序逐桶相等，四阶对照第1轮合并秩与距离两个完整边缘。

## 保护峰分母

```sh
python3 -B n10-protected/verify.py | cmp - n10-protected/verification.json
```

包含远端20、本机16个实际模数作业的全部原始模结果。用于正式生产的3个远端分隔是 5555555555、5445555566、4556555555，本机分隔是 4544555666、5445656555。每分隔并发至多4模，每项表数组预算512 MiB。已实际核验回执中的全部 PID，取回并复验后没有这批计算进程存留。重跑整数计数和前置小门槛的命令见 [四模计数器](n10-protected/COUNTER.md)。直切分母文件保留原生产回执的哈希，另附上整数分母数据，供可移植重放读取。

## 四象限完全交错

从 n10-gluing 子目录运行。两份 C++ 源码以同样的 clang 参数分别编译。每份程序分别用输入 2、3、5 生成相应秩表，六份新秩表应与归档 JSON 数据一致。两份源码均已实际重编译并比对通过。

```sh
python3 -B quadrants.py ../n10-protected/straight-upper.json \
  | cmp - n10-quadrants.json
python3 -B merge_quadrants.py "$SCRATCH/gluing-summary.json" \
  n10-quadrants.json ../n10-protected/straight-upper.json
cmp "$SCRATCH/gluing-summary.json" gluing-summary.json
python3 -B add_dominant.py "$SCRATCH/gluing-summary.json" \
  "$SCRATCH/dominant-summary.json"
cmp "$SCRATCH/dominant-summary.json" dominant-summary.json
```

合并器对小象限的完整排列和短相对词的交错进行交叉核验。大象限连接点秩表以及各距离门槛的原条件有理下界均在账本中。MacBook \(n=32\) Metal 批次独立运行，当前轮次只记录已观察到的进度。
