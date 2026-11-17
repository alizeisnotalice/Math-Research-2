# 四批有限 plug-in 一致性复核

本包以 Decimal(100) 后转 Fraction 的精确符号比较，分别检查四份已完成终端批的
\(B/A_2\le\sum_d dP_2(d)/[2(d-1)]\) 必要条件。为避免重复归档 r024/r025 raw，
审计器默认从同一 `certificates/` 根下的既有 r024/r025 证书读相对输入；必要时可用
`A022_CERTIFICATES_ROOT` 指定该根。n=16 MacBook 仍由既有 gzip 透明解压并核对原 raw SHA。

四个 slack 都为正，故这只说明它们未被该一个必要条件排除。它不证明原条件律、
混合、方差、显著性、代码正确性，也不允许跨批合并。
