# 复现

仅需Python3标准库，从本目录运行：

```sh
python3 -B replay.py --out /Volumes/data/CodexTemp/T2/a022-r036-replay-new
```

输出目录必须新建。脚本复制最小输入后重放8项小规模有理/枚举检查；不重复远程搜索、
n11图扫描或GPU生产。候选搜索的冻结源码与可选生成命令见finite-candidates/README.md。
工作文件名r042–r047是临时子任务号，公开轮号为36。
