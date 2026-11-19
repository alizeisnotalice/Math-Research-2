# 固定弧数 Reuleaux 多边形的半径亏损和稳定性

设 \(P\) 是有 \(m\ge3\) 条正长圆弧的宽一 Reuleaux 多边形，\(T_m\) 是同弧数的正形，\(m\) 为奇数。以 \(R,r\) 表示 \(P\) 的外接与内切圆半径，置 \(R_m=1/[2\cos(\pi/(2m))]\)，\(r_m=1-R_m\)。已有的定理给出 \(R\ge R_m\)，等号只在正形成立。宽一体的内外圆同心，\(r+R=1\)，故 \(P\) 到直径一圆盘在平移下的最佳 Hausdorff 距离恰为 \(R-1/2\)。以下直接给出角度证明，并使接近等号时的估计显式化。[Chen、Chen 与 Jin](https://link.springer.com/article/10.1007/s00022-023-00667-5) · [Henrot 与 Lucardesi](https://jep.centre-mersenne.org/articles/10.5802/jep.150/)

将外心取原点，按直径图的星形环顺序记角点 \(v_0,\ldots,v_{m-1}\)。每两个环邻点相距一，整个体为以这些点为圆心的单位圆盘之交。反过来，直径至多一的平面点集若有一条经过全部顶点的奇直径环，对应圆盘交是唯一的宽一完备体。[Kupitz、Martini 与 Wegner，定理 0.1 及 1.7](https://emis.de/ft/11518)

当 \(R<1\) 时各顶点半径在 \([1-R,R]\)。将相邻极角差提升成 \(\pi+e_i\)，其中 \(-\pi<e_i<\pi\)。环的奇性迫使 \(\sum e_i\) 为 \(\pi\) 的奇数倍。写 \(\alpha=\arccos(1/(2R))\)。距离平方对两个顶点的半径分别为凸函数，在 \([0,R]^2\) 上的最大值为 \(\max(R^2,4R^2\cos^2(e_i/2))\)。因此 \(|e_i|\le2\alpha\)，有 \(2m\alpha\ge\pi\)，即得 \(R\ge R_m\)。等号需要同向的 \(e_i=\pm\pi/m\) 且各半径相等，正星形的圆盘交即 \(T_m\)。共心圆盘的 Hausdorff 结论也可直接由 \(\min_u h_P(u)=1-\max_u h_P(u)=1-R\) 看出。

设 \(R<1/\sqrt2\)，且
\[
 B=2(m-1)\left(\alpha-\frac{\pi}{2m}\right)<\frac{\pi}{m},
 \quad T=\frac{4R^2\cos^2((\pi/m-B)/2)-1}{1+2R\cos(2\alpha)} .
\]
则模刚体运动有
\[
 R-R_m\le \inf_g d_H(P,gT_m)
 \le\frac{R_m}{1-R_m}\,[T+2R_m\sin(B/4)].
\]
为证上端，因 \(2m\alpha<3\pi\)，反向重排后可令 \(\sum e_i=\pi\)。各 \(e_i-\pi/m\) 的正部总和至多 \(B\)，负部总长相同。部分和的极差不超过 \(B\)，选一个公共转角后各顶点的极角误差至多 \(B/2\)。将一个顶点的半径记作 \(R-t\)，保持之固定，另一个半径取端点最大化距离，可得
\[
 4R^2\cos^2(e_i/2)-1
 \ge t[2R(1+\cos e_i)-t]\ge t[1+2R\cos(2\alpha)].
\]
这里 \(t\le 2R-1\)。又 \(e_i\ge\pi/m-B>0\)，所以 \(t\le T\)。由 \(4\cos^2(\pi/(2m))(R+R_m)\ge1+2R\) 可知 \(T\ge R-R_m\)。每个顶点距正形对应点至多 \(\Delta=T+2R_m\sin(B/4)\)。

正形各圆心模长 \(R_m\)。若 \(x\in P\)，它距每个正形圆心至多 \(1+\Delta\)，则 \(x/(1+\Delta/r_m)\in T_m\)，这是三角不等式。反向对 \(\Delta<r_m\) 有 \((1-\Delta/r_m)T_m\subset P\)。当 \(\Delta\ge r_m\)，用 \(0\in P\) 直接控制 \(T_m\) 到 \(P\) 的距离。下端由外半径对 Hausdorff 距离是 1-Lipschitz 得到。固定 \(m\) 时上式在 \(R_m\) 附近为 \(O_m(R-R_m)\)。对 \(m\ge5\)，可固定正形直径环中的 \(v_0,v_3,\ldots\)，小幅沿单位圆移动 \(v_1\)，令 \(v_2\) 取距 \(v_1,v_3\) 均为一的连续圆交点。其余点对的严格余量保留，三固定非共线点使一般小变形与原形不全等。故线性阶可在 \(R\downarrow R_m\) 的非正形族上检验。
