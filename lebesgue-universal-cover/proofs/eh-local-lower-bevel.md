# 下端真耳的局部化、面积与排除

在十一段 EH 的最后一条线段上记上端铰点为 \(X\)，下端为 \(T\)。取 \(Y=T+(X-T)/25000\)，沿 \(B_3\) 的真单位圆弧逆时针转 \(\alpha=10^{-7}\) 弧度到 \(D=B_3+R_\alpha(T-B_3)\)，加入从 \(Y\) 指向 \(D\) 的保留左半平面。

在 \(Y,D\) 的两处左转严格；新直线将 \(X\) 放在保留侧、将 \(T\) 放在删除侧。\(D\) 对旧十一条线的最小余量大于 \(4.59116\times10^{-8}\)。旧边界为凸曲线，直线 \(YD\) 已在 \(Y,D\) 各横穿一次，因而删去的真耳恰由线段 \(YT\)、真圆弧 \(TD\) 及弦 \(DY\) 围成。整个圆弧距 \(T\) 不超过 \(\alpha\)，而 \(|Y-T|<3.366193\times10^{-7}\)。由圆盘凸性，真耳的全部点都在此半径内。\(p=.427,\delta=.00205\) 的 AH 裁线在包住真耳的矩形上的余量超过 .46216。

真耳的面积为
\[
 \tfrac12[T-Y,D-Y]+\tfrac12(\alpha-\sin\alpha)
   =7.7273752630828941760593851285627551\times10^{-15}
\]
（显示末位按 Arb 区间取舍，完整有理外包见证书）。这项删除破坏普适性。由[强制透镜定理](eh-mandatory-lower-lens.md)，取 \(p=T,\lambda=10^{-8}\) 的定宽一体在 \(P\) 中只有一幅几何像。它含有 \((1-\lambda)T\)，该点落在新切线外的法向距离大于 \(3.0953974\times10^{-8}\)。

~~~sh
python -B computation/certify_eh_local_lower_bevel.py computation \
    computation/eh-upper-chain-parameters.json \
    computation/ah-height-switch-certificate.json OUTPUT.json
~~~
