# 第49轮复现说明

从本目录运行。用 data 卷的空目录保存重编译和比较结果。实测解释器 Python 3.15.0，编译器 Apple clang 21.0.0，整数算法无随机容差。

## 解析门槛

```sh
python3 -B suffix-inequalities/audit.py | cmp - suffix-inequalities/verification.json
python3 -B suffix-inequalities/audit_reflection.py \
  | cmp - suffix-inequalities/reflection-verification.json
python3 -B geometry-obstructions/audit_suffix_clusters.py \
  | cmp - geometry-obstructions/audit_suffix_clusters.stdout
python3 -B geometry-obstructions/audit_reflected_seam.py \
  | cmp - geometry-obstructions/audit_reflected_seam.stdout
python3 -B geometry-obstructions/audit_reflected_bands.py \
  | cmp - geometry-obstructions/audit_reflected_bands.stdout
python3 -B low-branch/audit_vulnerable.py | cmp - low-branch/vulnerable-verification.json
```

后缀小图验证用联合峰整数容斥与独立后缀递推两个方向。反射束用固定种子选前缀，Python 整数计算实际完成数。三份几何反例的测试使用显式坐标。分支修复的充分条件和反例都在对应解析证书中，不取决于随机频率。

## 低支与首接触的五维整数表

以下 SCRATCH 指向 data 卷临时目录。

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Werror low-branch/low_branch.cpp \
  -o "$SCRATCH/contact_counter"
"$SCRATCH/contact_counter" 3 4 8 | cmp - low-branch/n3-contact.tsv
"$SCRATCH/contact_counter" 4 6 8 | cmp - low-branch/n4-contact.tsv
"$SCRATCH/contact_counter" 5 6 8 | cmp - low-branch/n5-contact.tsv
python3 -B low-branch/verify.py | cmp - low-branch/verification.json
```

参数依次为尺寸、线程数、低支大小上限。标准输出的列依次为合并秩、低支大小、首接触间隔、原峰距、其他外邻数、词数。标准错误给出每任务的资源计划。五阶正式五维表使用6个线程，数组预算上界3,372,220,416字节，实测峰值驻留内存2,999,779,328字节。四维投影另与独立12线程运行逐项核对。该五维表不包含大小超过8的低支。

## 十一阶整数分子和两种分母

```sh
python3 -B odd-quadrants/verify.py | cmp - odd-quadrants/n11-numerator.json
python3 -B odd-quadrants/protected_upper.py | cmp - odd-quadrants/four-block-bound.json
python3 -B n11-dense/verify.py | cmp - n11-dense/verification.json
```

奇数四象限的矩形单峰根数重新做了三模恢复及独立 Python 全整数核对。重算三个区域模数的办法见第46轮 layered/REPLAY.md，轮廓取 555555，表预算每项64 MiB，合计192 MiB。另两张方形根账引用本仓库整数文件。

稠密五模分母计数器见 [COUNTER.md](n11-dense/COUNTER.md)。实际生产按两个区域分侧。把 n11-dense 的三个源输入文件与 run_side.py 复制到各自空目录，采用 COUNTER.md 的 Apple clang 命令编译。对本机轮廓 55555566666 与 MacBook 轮廓 66666655555 各按索引0、1、2、3、4依次运行。

```sh
# PROFILE 是该侧轮廓，BUDGET 本机4480、MacBook6656
python3 -B run_side.py "$PROFILE" "p$INDEX" --max-table-mib "$BUDGET" \
  --prime-index "$INDEX" --timeout 10800
```

归档 complete-batch.json 保存十个原始模数标准输出、资源回执和驱动终态。收回后对全部原始文件做哈希验证，核对记录的待回收 PID 均消失。verify.py 重验5个素数、完整模积、独立集层数、CRT、每侧阶乘边缘及全图两类最高秩扣项，并用同一严格分子给出原条件律精确有理下界。有限证明如第49轮报告所述。
