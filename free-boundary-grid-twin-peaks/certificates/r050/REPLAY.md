# 第50轮复现说明

使用 data 卷的空目录作 SCRATCH。解析及整数核验使用 Python 3.15.0 和 Apple clang 21.0.0。以下命令从本证书目录运行。浮点诊断另有 Swift 版本、种子域和实际回执。

## 首次合并低支删除

```sh
python3 -B delete-geometry/audit_deleted_suffix.py \
  | cmp - delete-geometry/audit_deleted_suffix.stdout
clang++ -std=c++20 -O3 -Wall -Wextra -Werror delete-branch/delete_branch.cpp \
  -o "$SCRATCH/delete_branch"
"$SCRATCH/delete_branch" 3 4 8 | cmp - delete-branch/n3.tsv
"$SCRATCH/delete_branch" 4 12 8 | cmp - delete-branch/n4-b8.tsv
"$SCRATCH/delete_branch" 5 12 8 | cmp - delete-branch/n5-b8.tsv
python3 -B delete-branch/verify.py | cmp - delete-branch/verification.json
```

三个位置参数依次为方格边长、线程数、低支大小上限。标准输出依次为首次合并秩、低支大小、原图峰距、全部词数、可删除词数。各低支按正方形刚性对称取代表。固定低支后，原图与删除它的诱导图分别建立合法增长后缀表，两侧不相混用。根距频数和合法穿插乘到两个表。输出脚本用旧的原条件律五维首接触账本验证全词边缘，三阶另直接从洞图单峰词枚举回插，验证全部可删除计数。五阶资源见 n5-resource.txt。

## 十一阶前缀条件计数

源码及公式在 [n11-gluing/PREFIX_GLUE.md](n11-gluing/PREFIX_GLUE.md)。使用该说明重编译，阶段性审查如下。

- 复制 n11-gluing 的程序及两个小计划、主计划到空目录，按说明编译和运行。将环境变量 A022_RESEARCH_ROOT 指向仓库中的课题根目录，以读取第48、49轮已冻结的输入。若就地读取归档数据，可省略该变量。
- 复编译后已重新求值全部4,150个主要上下文，与归档逐项一致。资源回执时间和 PID 可以变化。判断一致性时比较各查询的参数、完成数、独立集层状态数和总查询量。
- 用主计划及原完成账重新运行 merge_n11.py，再运行 add_dominant.py。两份整数分子、重叠和条件有理下界逐字节重得。分母直接引用第49轮的完整五模保护峰上界。两个小门槛分别验证 \(3\times3\) 的 \(9!\) 及五阶450个秩桶。

## 七阶 Metal 有限诊断

```sh
python3 -B local-metal/verify.py | cmp - local-metal/summary.json
```

在空目录中以 Swift 优化模式编译第26轮校准源，使用共享 data 模块缓存，实际参数依次为 7、65536、16、2、2、输出路径、9e50128000000001、strat、128。JSONL包含完整16个重复、完整秩提升汇总。校验器比对原始文件、资源和源码哈希，核对逐重复种子、整数矩关系与进程终态，再用本仓库的原律整数真值校准。诊断的有限接近程度不提供收敛保证。
