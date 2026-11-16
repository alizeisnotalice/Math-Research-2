# 六点总次数44：完整精确账本

新增43/44层共43个有序度数对、22片：912953731个原始支撑对，441939469个规范代表；439278131个有单位碰撞证书，2661338个唯一表示公平，未知0。结合此前≤42及常数因子/共同指数约化，严格范围为每侧≤6点、总次数≤44。

每片 `verification.json` 来自独立逐项验证器，`job.json` 确认进程退出。主线程 `check_six44_receipts.py` 又解压全部账本，核对实际内容哈希、完整度数覆盖、精确组合计数、源代码哈希和无余例；27.694秒通过。完整可再生gzip账本留在CodexTemp，不复制数十GB至Git。

复现每片：用课题 computation/support_search_layer.py 按summary.config生成，随后 `python3 computation/verify_support_layer.py <该片目录>`。聚合：`python3 check_six44_receipts.py --ledger-root /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round024-fullcpu/output --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/six44-recheck.json`。本聚合不替代已完成的逐项重算。
