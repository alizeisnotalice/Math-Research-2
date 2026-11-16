# d=11 前缀池与 CPU-only 准备

独立标准库枚举了全体 `3^10=59,049` 个低前缀编码，不导入或运行 GPU
worker。结果为 5,131 个合法低前缀、3,181 个通过两项已声明必要过滤；
SHA-256 分别为 `db4403d00b7fdd18b13a30ac7e7f76f0191f86e213b95a0549ef62eb8f6a1021`
和 `821c24a85b9d6597cf5714dfab0ecd3cc9b6775da4b8f53e070903d2ca25d15f`。

`worker/` 是 round028 d10 frozen 三个源的独立副本。仅把 worker 的两处
`d<=10` 门改为 `d<=11`、错误文本同步改为 `1..11`，并把 driver 固定的
worker SHA 更新为新值；详情见 `source.diff`。运行 CPU-only 的
`m=48,n=88,d=11,starts=32768,steps=12000,restart=2000,margin=.05,seed=2026111609`
结果与独立池逐项一致，且 starts 覆盖 3,181 个前缀。

`margin=.05` 是有意只探索 `t∈[.05,.95]` 的数值子模型，不能排除贴近
`t=0,1` 的情形，更不能作为一般排除定理或精确反例证书。没有启动 GPU/MPS。
