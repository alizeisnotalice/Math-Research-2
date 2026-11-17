# 来源及核验范围

- 用户附件：A022_Pro_原问题闭合提示词_20261115.md；来源路径 /Users/zzh/Downloads/A022_Pro_原问题闭合提示词_20261115.md；SHA-256 `c7b46a6422bff416afe5ee06c54e7c80ff0b29b08b7fb532b58851a60110b0bc`。已阅读全文；“已证”描述不自动作为本轮证明。只归档自写推导及来源索引，不复制附件全文。
- Burdzy, Problem 9：https://sites.math.washington.edu/~burdzy/open_mathjax.php 。2026-11-15 实查，第9题给出 n×n 方格、两峰条件和 ρ_n/n→0 问题；按本任务固定为自由边界。
- Burdzy–Pal, Twin peaks, arXiv:1606.08025v2：https://arxiv.org/html/1606.08025v2 ；期刊 RSA 56 (2020), 432–460。已核对摘要、§1、§2、§3定义：§1 的否定猜想为环面；§2 Example 2.1 说明均匀边界生长不等于条件均匀标号。第14轮于2026-11-16进一步实际核读§4 Theorem 4.2及完整证明（HTML行99–149），其定理对象是自由方格，使用其中固定根单峰数下界；双峰转移由本项目另行推导，见proofs/two-peak-level-boundary.md。不归档论文全文。
- GitHub 文件限制：https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github 。2026-11-15 实查，50 MiB 警告、100 MiB 阻止；项目执行更严格的单文件 10 MiB 默认。

- 第2轮独立子任务初始面积阶界：/Volumes/data/CodexTemp/T2/a022-r-lower/r_lower_bound.md；SHA-256 `95ac5006706fd4b48c45faabf3ec5afff893629f1d8d8853072fb03f7ba359f3`。已逐步核验提升合法性、层析计数和逆像；主Agent进一步推导逐输出权重、远峰下界及精确逆像公式，见proofs/weighted-promotion.md。

- 第18轮于2026-11-16重新核读 Twin peaks v2 §1–2、§4–5；继续区分环面猜想、自由方格粗糙边界定理和固定宽矩形。新根矩/可拆分支/前缀核结论为本项目自写推导，不声称来自该论文。


第42轮网页复核（2026-11-17）：重新读取 [Burdzy Problem 9](https://sites.math.washington.edu/~burdzy/open_mathjax.php) 的原问题表述和 [Twin peaks arXiv摘要页](https://arxiv.org/abs/1606.08025)。对象仍为均匀标号在恰两峰下的n×n方格图距；本轮仅重新核对网页与摘要，不把摘要中的树图结论当方格结论，也不声称本轮重审了全文证明。构造和有限证书是本项目自写。
