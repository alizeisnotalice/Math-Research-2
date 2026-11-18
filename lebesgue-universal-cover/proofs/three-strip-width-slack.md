# 三夹带两帽阻挡的宽度余量加强

以精确三夹带的两谐波阻挡为基线。保持其中帽角、帽深和角扰动的同一范围，唯独允许每条居中夹带的半宽从 \(1/2\) 增加到 \(1/2+\eta_i\)，其中
\[
 0\le\eta_i\le10^{-5}M,\qquad
 M=\max(|\epsilon_1|,|\epsilon_2|)>0 .
\]
既有两个光滑定宽体中至少一个仍无法放入。因而这组帽法向与深度内的普适覆盖，至少需要一条夹带的全宽严格大于 \(1+2\cdot10^{-5}M\)。

设某一候选等距像在 \(n_i\) 上的有符号中心误差为 \(s_i\)。定宽一使 \(|s_i|\le\eta_i\)，而归中平衡变成
\[
 \rho_k\left|\operatorname{Re}(e^{-ik\phi}S_k)\right|
 \le\sum_i\lambda_i\eta_i\le3\cdot10^{-5}M.
\]
按基线相位估计，
\[
 |S_k|\ge (k-1)(3M/4)(1-r_k),\qquad
 r_k=(k^2+2k-1)10^{-3}/[(k-1)(3/4)].
\]
因此相位 \(e^{-ik\phi}\) 与两种精确相位 \(\pm\overline{V_k}/|V_k|\) 之一相差少于 \(0.5^\circ\)，有理证书核验了对应的 \(\sin(0.5^\circ)\) 正余量。

以基线同一方式从 \(k=7,13\) 中选择一个 k。额外相位误差扣除后，两个帽残差的谐波投影仍异号，其绝对值均大于 \(|Z_k|\sin1.5^\circ\)。夹带余量导致帽支持增量还多了 \(a_\beta s_0+b_\beta s_1\)。因 \(|b_\beta|\le2\)、\(|a_\beta|\le3\) 且 \(M\le10^{-3}\)，绝对误差至多 \(5\times10^{-8}\)。扣除后仍有一个帽被击穿超过 \(1.1055\times10^{-5}\)，严格大于允许的 \(10^{-5}\) 外放量。

独立入口 computation/certify_three_strip_width_slack.py 实际以 Python 3.13.18 执行退出零。它重放原精确夹带的有理三角检查，再核验本加强的额外相位和支持余量。computation/three-strip-width-slack-certificate.json 记录两份源码哈希和精确余量。
