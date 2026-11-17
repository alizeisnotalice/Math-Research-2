# 第 53 轮 \(Q_3\) 公式文字更正

更正前文件 `s234_low_reduction_zh.md` 的 SHA-256 为
`defc580eeaf33d9029c1f974abd6c3ddc01cc44d5d8864969d5ca87a08e010b8`。
第 36 行曾误写为 \(G_3=Q_3=-\ell G_2-k\)。

正确的双侧第 3 行是
\[
 G_3=-\ell G_2-k,\qquad Q_3=-\ell Q_2-k.
\]
除非 \(G_2=Q_2\)，两者不相等。错误只在未发布说明文字和向主线程的摘要中。

冻结的 `s234_core_box.py` 只使用 \(G_3=-\ell G_2-k\) 作记录字段，且其
\(s_2,s_3,s_4\) 核心计数不使用任何 \(Q_3\)；故已运行的初始盒输出不受影响。
独立双侧形式卷积脚本 `s234_bilateral_identity_check.py` 另外验证了在
\(G_1=Q_1=-\ell\)、\(b=\ell a+k\) 下，\(A_3=B_3=0\)。
