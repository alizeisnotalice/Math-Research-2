# d=12 前缀池与 CPU-only 准备

独立标准库枚举全体 `3^11=177,147` 个低前缀编码，得到 11,005 个低端合法编码和
6,605 个通过两项已声明必要过滤的编码；SHA-256 分别为
`0499ace4436dcb94dd12930609f086992aba76fb45fe4d5291f5ea3d8346a5a3` 和
`b40d0e2f03d5c92b50c832d695c874fd2b77fe416c1ca2af846c829b1c245ad8`。
计划 `starts=32,768` 大于 6,605，故覆盖算术门通过。

`worker/` 是 round030 frozen d11 的独立副本；仅将 worker 的两处 `d<=11` 改为
`d<=12`、错误文本改为 `1..12`，并更新 driver 固定 worker SHA。实际 CPU-only
运行参数为 `m=52,n=92,d=12,starts=32768,steps=12000,restart=2000,margin=.001,seed=2026111610`，
其 pool 结果逐项匹配独立枚举。没有启动 GPU/MPS；该有限输入准备不构成任何排除、
存在或精确反例结论。

主线程另从归档副本实际复跑独立枚举（session85690退出0）与 worker CPU-only（session9618退出0），全部计数和两份池哈希吻合。`parent-*.json` 保存实际回执。冻结 worker 中的 `FROZEN_SOURCE` 绝对路径仅为可选来源哈希元数据：不存在时记 null，程序不导入或执行它。运行依赖为同目录 `fractional_prefixes.py`、NumPy、PyTorch及可用MPS。当前 d11 尚在运行，d12 未启动。
