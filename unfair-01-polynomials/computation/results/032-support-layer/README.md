# 42／34层完成回执与精确余例

42／34层全部38个搜索/独立验证分片结束，session39103退出0、父3440已退出。六点部分21片覆盖总次41/42的41个度数对，无余例，因此与已有≤40范围合成 [六点总次≤42定理](../../../proofs/032-six-point-total42.md)。九点部分仍有1个未解决余例，严格范围保持≤32。

`six42-m*/`保存每片参数、计数、独立验证回执及空的余例文件。完整逐项账本约15.2 GB未压缩内容可再生，保存在 `parent-six42-coverage.json` 所列 data 临时路径；未把大量可再生原始账本重复写入Git。该索引记录每份压缩与未压缩SHA以及字节数。原始文件仍可逐项读取，或使用已提交的搜索器重新生成。

本轮主线程完整解压并核对21份账本哈希、源码绑定、全部次数对及组合原始计数，运行16.99秒。逐项数学验证由各分片的独立 `computation/verify_support_layer.py` 已执行并记录；不把本轮哈希检查说成再次执行全部数学判定。

```sh
python3 check_six42_receipts.py --ledger-root /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round020-fullcpu/output --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/six42-integrity-recheck.json
python3 verify_residual_certificates.py
```

若需从头重跑某片，用其summary中的config调用 `../../support_search_layer.py` 输出到新临时目录，再调用 `../../verify_support_layer.py <该目录>`。各片采用整数算术，不涉及实系数的浮点采样。

`residual-46-mapping.json`保留九点层全部46个余例及39个严格对角支配见证。`parent-residual-verification.json`由主线程从支撑重新构造对角项并验证注入、严格性，未调用搜索器寻找见证。其余7条中又有6条取得正多项式恒等式，主线程执行 `verify_positive_six.py`，用不导入搜索器的Fraction验证器重建完整卷积方程、严格变量范围及正恒等式。45条现已认证，余1条没有证书，不能把匹配失败当排除，也不能称作原题反例。运行中的44／36层不属于本包的完成范围。
