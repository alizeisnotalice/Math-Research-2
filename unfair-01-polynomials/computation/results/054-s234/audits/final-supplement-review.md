# 第54轮最终补充审查：正式稿、full 回执与身份链

审查日期：2026-11-17（Asia/Shanghai）。本补充覆盖此前审查完成后写入的 full
回执、正式稿文字修订和独立降幂根筛；不重跑 499 秒封包。

## 被审对象与实际身份

- 正式证明：`sync/unfair-01-polynomials/proofs/054-degree7-complete-quadratic-core.md`，
  SHA-256 `77df9abd64f0ddae249a0f8b1f16884b535ed34e67ed60e58aaf975feaaea5b0`。
- full 回执：`round054-parent/packaged-full-final.json`，SHA-256
  `67bff5e9c909c905e4c233bc2b3541602e3a0df2869590a85a59924d502937b7`。
- 核心清单：`round054-parent/package/core-manifest.json`，SHA-256
  `64d5aeeaa701ad969253ddf797f4779da7e2dbe11a3e362e748ee7478306271c`。
- 独立降幂根筛总回执：`round054-parent/descending-full/summary.json`，SHA-256
  `b5a27ad9ace02160f26450444c10931791d7a291e40f3031ef4225a53031931f`。

正式稿与最终 parent draft 在本次检查时逐字节相同。此前审查中“七个正轴点”是
历史文字笔误；正式稿第78行和实际独立筛均使用 **九个** 点：

\[
1,\tfrac12,2,\tfrac14,4,\tfrac34,\tfrac43,\tfrac65,\tfrac56.
\]

## 最终数学文字复核

1. 正式稿第27行现在明确取倒数根 \(w_i=1/z_i\)，写
   \(P(x)=\prod_i(1-w_i x)\) 与 \(s_j=\sum_iw_i^j\)。因此第31--34行是按
   x 升幂写出的 Newton 恒等式，方向正确；倒数根仍在同一严格环带且模积为一，
   所以 \(T_j\) 界没有改变。
2. 第62行已收紧为“允许自由整数 \(G_6,Q_6\) 平移后”的说法，正确区分了
   `rho,sigma` 的整数据依赖和实际单位区间系数的无理部分。
3. 第114行已改为“G、Q、P 均首一且常数项一”，消除了原 `P/G/Q` 的歧义。
4. 全文仍明确：根/正轴/Schur 均为单向必要筛；最终无界结论来自完整空树与有限
   余因子补零，而不是把余因子次数截到71。常数 `G/Q`、带符号中间系数、以及
   `P,Pbar` 不互素都在量词内。

未发现本次修订引入的数学漏洞。

## 回执与绑定审计

我实际运行独立只读脚本 `verify_final_identity.py`（SHA-256
`250a81cfb0184153515ae3a72510867a75b9a01a7dc10d73a299e09bc0e26273`），实际回执为
`python-math-final-v3/20261117T141026Z-82ac27c9/run.json`。其输出
`final-identity-audit.json`（SHA-256
`e0128efd37f56ed97bc1ad77c84e7aa2e57827e8eb121dc86bf0d8923e74552b`）完成：

- 重新计算核心清单所列 **212** 个文件的字节数与 SHA-256，全部匹配；
- full 回执的 core-manifest SHA 与实际清单相同，driver SHA
  `b2ad7be4f5874dce2f3a84e986f7c37f4f95203cfe92ca293b93d370a4c0b1cd`
  与封包中的 `verify_package.py` 相同；
- full 实跑状态为 passed，`fields=58`、`cached_fields=0`、`direct_AB_replayed=true`、
  `54,726,220` 个母式根码、`68,593` 棵空树和 `621,704` 树层；耗时
  `499.828648750` 秒；
- 独立 descending 根筛的58份逐域 JSON 均为 passed，并按实际读取的
  `round053-parent/row567-full` 与 `round053-parent/annulus-full` 输入/父账本逐份重算
  SHA；逐域合计与 descending summary 及 full receipt 的根筛统计完全相同；
- descending 回执中的源 SHA 对应冻结的 original 源码
  `verify_s234_annulus_descending.py`。封包内 adapted 版本与 original 的唯一 diff
  是 CLI 默认 I/O 路径，非根筛算术；这一点也与 `io-adaptations.json` 一致。

这是一条“实际 full 重算 + 独立逐域降幂根筛 + 本次清单/绑定重哈希”的证据链，
不是仅凭哈希相等登记计算成功。

## 局部 r=1,t=5 补充

未扩大局部搜索。另一个独立 Fraction/完整整数环验证器
`round054-parent/local-r1t5/verify_fraction.py`（SHA-256
`2fbbc49bc044b063f24e2c0dcd3ecb44005aa7780062195a965a6de31bbffb3d`）的
`fraction-audit.json`（SHA-256
`1df06be876f0eeb72f02af62c9236dd95733d8c7213ad2c5603aacc0d2dbdd0c`）记录：所有
局部非零非单位系数为代数整数，\(C_0,\ldots,C_5=1\)，而
\(C_6\in(0,1)\)。其 scope 明确只是局部前缀，绝不是完整分解或原题反例。

## 审查结论

正式稿现可将“七次实二次原始块的 S={2,3,4}（及反转 S={3,4,5}）”作为已完成的
计算辅助排除定理登记；结合已审第43--52轮的边界，非有理实二次原始块必须有
\(\deg P\ge8\)。该结论不推广为一般均匀和猜想，也不覆盖非二次系数域。
