# 计算与复现

在本目录运行以下命令可复核五维证明中的有限表及六变量扩张结果。参数 Walsh 表使用整数计算，投影界的比值用有理数比较。

```sh
python3 parameter_certificate.py
python3 audit_five_minimum.py
python3 check_second_weight.py parameter_certificate.json second-weight-audit.json
python3 verify_tower.py
```

六变量扩张枚举使用 `cc -O3 tower.c -o tower`。命令行依次传入分片号（0至107）、108、第38轮第六个字的十六进制掩码 `000000003c5a3c5a`、此前选定的所有更高维字。下表各行是沿指定路线贪心选出的下一字，与其支持并的最小点数。每次均完整覆盖六变量三次及以下多项式。不同路线的极值可能不同。

|新维数|新字掩码|支撑|可选独立方向数|
|---:|---|---:|---:|
|7|000000000f0f0f0f|32|229312|
|8|0066006666006600|40|229248|
|9|00aa00aa00aa00aa|44|98048|
|10|0033003300330033|48|32256|

本机12个进程和 MacBook 15个进程分别处理余数0至47和48至107。Python 独立验证了各级每个非空分片的见证。九维批次的 SSH 连接一度超时，保留本机已完成分片，重连补齐远端，最终每级均有108个互异分片。

第38轮的CPU分片另由两台机器的 Apple Metal 检验。编译 `swiftc -O WitnessMetal.swift -o witnessMetal`，分别以 `../round038/mini-complete.txt` 和 `../round038/macbook-complete.txt` 为参数。每个非空分片取其见证，GPU 对64个相位逐个计算64项 Walsh 变换。本机 Apple M6 执行24个命令缓冲区、核验98304个系数；MacBook Apple M5 Pro 执行30个、核验122880个系数，全部绝对值为8。Metal 时间分别约0.0011秒和0.0416秒。GPU短批的完成能由命令执行状态确认。

11月18日09:04（上海），第38轮CPU分片同步运行时，本机两次测得整机CPU空闲0%，MacBook五次；A28这批进程合计峰值采样分别为1081.3%和956.5%（单核100%）。同期GPU设备利用率采样100%为整机数据，多课题共享。
