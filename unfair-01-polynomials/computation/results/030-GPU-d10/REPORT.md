# 第 30 轮：d=10 完成 GPU 批次的独立审计

审计对象是 `round028-prefix-gpu/live-m44n84d10`，参数为
`m=44,n=84,d=10,starts=32768,steps=12000,restart_every=2000,margin=.001,seed=2026111608`。
保存 config 的实际软件记录是 PyTorch `2.15.0`、NumPy `2.5.4`、MPS、float32、
macOS `27.2-arm64-arm-64bit-Mach-O`。Worker、driver、prefix module 的 SHA-256
分别为 `b7b410d1339c3e72ca1a01d73e00fbf3dc922484e5226522ca1d8f2a39ccced0`、
`508d79530efb8a4b252a35704decd7afd0deed298ec69d6cb0916461f32575ce`、
`0b1a06770a61fe113031209986f676c7e1afb862d27c0c697dbebc0a9f78481d`，均与 config/
queue 记录匹配。

独立标准库覆盖账本重新枚举 `3^9=19,683` 个 d=10 前缀：2,463 个低端合法、
1,671 个通过两项已声明必要过滤（SHA 分别为
`c77752411ba366a58dce63c939b3aa640b932f2661d4ce77cd689402978ff6fe` 和
`a93e50782e8a105c4f1d452d6d96c0b533960ef6e2c78829661a435446892641`）。冻结 worker
的确定性轮换赋值式经源 SHA 和文本复核。六个 complete cohort 均有 32,768 rows：
每个前缀在每 cohort 出现 19 或 20 次（652 个为 19、1,019 个为 20），所以每 cohort
覆盖全池；计划总行数 196,608。五份 restart（steps 2000--10000）与最终 step 12000
checkpoint 都存在，故这是完成账本而非 partial。

独立 verifier 不导入 worker/prefix module，以 `Fraction` 对全部 5×16+16=96 个保存
JSON 样本重卷积，并检查端点、非负、margin、assigned/derived prefix、低卷积及 d-cross。
所有行通过；96 行互异，partial rows 为 0。最小保存的精确 Fraction residual 是
`0.0020387854892760515`，与保存 float residual 的最大差为
`2.988241283390671e-19`。因此 **没有任何保存样本是精确二元乘积/反例**。
这不声称未保存的浮点轨迹被精确认证，也不把有限前缀或数值阴性提升为一般排除定理。

`queue.json` 为 complete、returncode 0 且 child reaped；`round030-parent/d10-exit.json`
记录 session 22180 exit 0、parent 27324/child 27325 已消失。退出前 release 元数据记录
97,256,192 bytes，仅描述退出前值，不能解释为进程结束后的常驻内存。

可复现材料：`audit_d10_coverage.py`、`d10-coverage-receipt.json`、
`verify_prefix_fraction.py` 和 `d10-prefix-fraction-receipt.json`。这次审计只读已完成
批次，未启动 GPU，未触碰运行中的 d11。
