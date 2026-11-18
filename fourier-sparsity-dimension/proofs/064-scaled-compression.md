# 仿射限制下的核压缩与纪录谱多尺度约束

设 \(F:X\to\{\pm1\}\) 有谱支撑 \(S\)，\(s=|S|\)，\(d=\dim\langle S\rangle\)。任选频率子空间 \(H\le X^*\)，先删去无关方向，可设 \(H\le\langle S\rangle\)；一般情形改用 \(h=\dim(H\cap\langle S\rangle)\)。以下记 \(h=\dim H\)，\(Y^*=X^*/H\)，\(n_\beta=|S\cap\pi^{-1}(\beta)|\)，\(r=\dim\langle\pi(S)\rangle=d-h\)。同一商标签中的频率称为一带。

**压缩定理。** 每个 \(M=2^k\ge1\) 都有一个布尔仿射限制 \(F'\) 满足
\[
 s'\le C_M=\sum_\beta\min(n_\beta,M),\qquad
 d'\ge r-\min_{I}\sum_{\beta\in I} e_\beta,
 \quad
 e_\beta=\min\Big\{\frac{n_\beta-1}{2M},\
              \frac{n_\beta-1}{M+n_\beta-1}\Big\}.
\tag{1}
\]
\(I\) 在所有由出现的商标签构成的线性基间取值。若实际知道该带的物理活动率 \(p_\beta\)，可将第二项替换成 \((1-p_\beta)/(1+(M-1)p_\beta)\)。

选择直和坐标 \(X=H^*\times Y\)。在 \(z\in\mathbb F_2^k\) 上以均匀的随机仿射映射 \(\phi(z)=Tz+b\) 替代 \(H^*\) 中的物理变量，保持 \(y\in Y\) 不动，得到 \(F'(z,y)=F(\phi(z),y)\)。允许 \(T\) 秩亏，去掉无关 \(z\) 即为通常的仿射限制。频率 \((a,\beta)\) 变为 \((T^*a,\beta)\)，系数附符号 \((-1)^{a(b)}\)，故每带至多 \(\min(n_\beta,M)\) 点。商标签严格保持分离。

对一个选中的 \(\beta\)，固定带内一个原频率 \(a\)。每个其余原频率与其同桶概率为 \(1/M\)。给定 \(T\)，若桶内有 \(j+1\) 个原系数，随机平移 \(b\) 下该桶非零的概率至少 \(1/(j+1)\)。这是基本不确定性界：一个有 \(u\) 个非零归一 Walsh 系数的函数，其物理活动率至少 \(1/u\)，由反演、Cauchy 与 Parseval 直接得到。故整带丢失概率不超过 \(\mathbb E[j/(j+1)]\le\mathbb E[j]/2=(n_\beta-1)/(2M)\)。

另令 \(g_\beta(x)=\widehat{F(x,\cdot)}(\beta)\)，其物理活动率 \(p_\beta\ge1/n_\beta\)。整带消失等价于 \(g_\beta(\phi(z))\) 在全部 \(z\) 上为零。不同 \(z\) 的 \(\phi(z)\) 两两独立均匀（\(M=1\) 时只用一个均匀点）。计数器 \(Z=\sum_z1_{\{g_\beta(\phi(z))\ne0\}}\) 有 \(\mathbb E Z=Mp_\beta\)，\(\mathbb E Z^2=Mp_\beta+M(M-1)p_\beta^2\)。Cauchy 给 \(\Pr(Z>0)\ge(\mathbb E Z)^2/\mathbb E Z^2\)，得到 (1) 的第二项。幸存的基标签始终线性独立，故对最优 \(I\) 在同一次限制上求期望，再取至少达到期望的一次实现即可。无需各带丢失事件独立。

**单平坦版本。** 对仿射核心 \(L\) 用其方向 \(H\)，写 \(q=|S\setminus L|\)。其压缩后至多 \(M\) 点，故 \(s'\le q+M\)。若 \(L\) 是线性空间，基中不会出现零标签，故尾部各带都可使用，且 \(\sum_{\beta\in I}(n_\beta-1)\le q-r\)。若 \(L\) 的公共商标签非零，可在商掉其一条直线之后选择代表，并避开该标签本身，得到大小至少 \(r-1\) 的独立集 \(I\)，有 \(\sum_I(n_\beta-1)\le q-|I|\)。因而无须对核心系数作幅度和相位假设，即有
\[
d'\ge d-h-1-\frac{q}{2M}.
\tag{2}
\]
若直接采用尾标签线性秩 \(r_{\rm out}\)，可改为 \(d'\ge r_{\rm out}-(q-r_{\rm out})/(2M)\)，不会因核心为仿射而额外减一。如果尾标签本身已经张成全部商标签，式 (2) 中的 1 可删去。

**纪录谱的截断质量。** 称 \(F\) 为纪录谱，意指所有 \(s_G\le s\) 的布尔函数都满足 \(d_G/\sqrt{s_G}\le d/\sqrt s\)。一般 A28 若失败，可沿纪录子序列令此比值趋于无穷。记 \(\Delta_M=\sum_\beta(n_\beta-M)_+=s-C_M\)。由 (1) 对每个 \(H,M\) 有
\[
1-\sqrt{1-\Delta_M/s}
 \le \frac{h+\min_I\sum_{\beta\in I}e_\beta}{d},\qquad
\frac{\Delta_M}{s}\le\frac{2h}{d}+
        \frac{\min_I\sum_{\beta\in I}(n_\beta-1)}{Md}
 \le\frac{2h}{d}+\frac{s-r}{Md}.
\tag{3}
\]
若右端的扣减小于一，左式还可平方保留二次修正。在所有 \(h=o(d)\) 的方向上，只要 \(M\gg s/d\)，全部陪集超出阈值 \(M\) 的总超额都是 \(o(s)\)，同时不要求核中本来有一个完整平坦。特别当 \(M\ge T s/d\) 时，均匀谱点 \(\alpha\in S\) 所在带的规模大于 \(2M\) 的概率至多 \(4h/d+2/T\)。

单一任意 \(h\)-维仿射平坦中有 \(s-q\) 个谱点时，由 (2) 与纪录条件得
\[
\frac{s-q}{s}\le\frac M s+\frac{2(h+1)}d+\frac q{Md},
\]
在 \(q+M\ge s\) 时右边也显然足够。取 \(M\) 为不小于 \(\sqrt{sq/d}\) 的最小二次幂，使最后两项与首项中的压缩误差之和至多 \(3\sqrt{q/(sd)}\)（\(q=0\) 时直接 \(d\le h+1\)）。故当 \(d\to\infty\)，任何 \(h=o(d)\) 平坦所占比例趋于零。此处只对纪录代表断言反集中；第62轮构造的稠密核心族无需刷新全体函数的比值纪录。
