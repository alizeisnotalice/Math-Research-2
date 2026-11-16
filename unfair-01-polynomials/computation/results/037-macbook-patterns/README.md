# MacBook q23–30：必要模式与完整逐项验证

远端15核/24GiB，已有SSH主机密钥严格校验。Python3.12.13，发现时15进程；第二批采用逐行gzip以限制内存。全部计算在RAM，结果流回Mac mini的data盘。
q23–26共27参数任务24.752秒，独立15进程完整位掩码/事件核验10.512秒；q27–30共29参数任务266.889秒，另由6进程流式验证全部56组中的后29组。完成的父子进程均已退出并回收。

ledger-index.json给每个实际原始gzip的SHA、字节数、覆盖和未知数，all-independent-receipts.json给逐参数完整独立核验。不把unknown当作可行分解。本轮解析下降不依赖这些有限结果，且不将它们计入原猜想的支撑次数排除范围。
原始数百MiB账本可再生，保留在CodexTemp对应目录；公开提交参数、算法、逐项核验代码、哈希、回执与小型回归，不重复归档全部轨迹。传播程序不限制母因子任意次数，因其只使用已证明的低端范围；然而参数枚举本身只覆盖这里列出的q。

可从本目录在data临时目录重新生成一个参数并独立重放：

```sh
/Users/zzh/.venvs/main/bin/python generate_one.py --p 2 --q 23 --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay-p2q23.json.gz
/Users/zzh/.venvs/main/bin/python verify_pattern.py 2 23 < /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay-p2q23.json.gz > /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay-p2q23-receipt.jsonl
```

新生成文件使用统一流式格式；原首批不是rows优先的格式，不能直接交给流式解析器。其原核验使用first-verify_remote.py。运行时间/PID以及JSON字段顺序使复跑字节不必等同历史SHA，数学覆盖与事件逐项验证独立成立。生产流式验证器已通过7类篡改拒绝与正例；stream-negative-tests.json记录。
