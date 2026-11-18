# 分带系数固定的商标签提升与外部秩

保留前稿的 $F(y,x)=\sum_{\beta\in B}a_\beta(y)\chi_\beta(x)$。$a_\beta(y)$ 为实数；记外部谱 $S_\beta=\operatorname{supp}\widehat a_\beta$，全稀疏度 $s=\sum|S_\beta|$。以下在限制到每个物理切片之后使用第5轮零和分块法，所有切片须使用同一组提升后的商标签。

令 $E=\mathbb F_2^B$。对各 $y$ 和原非零差 $\delta$，将 $\{\beta,\gamma\}$ 满足 $\beta+\gamma=\delta$ 且 $a_\beta(y)a_\gamma(y)\ne0$ 的无序配对分成若干组，各组产品之和均为0。因切片为布尔，整批本身已是一个允许组。对同组中任意两对，取其端点指标的异或，在 $E$ 内张成 $R$。定义 $u_\beta=e_\beta+R$。原标签映射 $e_\beta\mapsto\beta$ 湮灭 $R$，故 $u_\beta$ 各不相同，并可线性投回原标签。

在新数据对偶空间 $E/R$ 中定义
$$
\widetilde F(y,z)=\sum_{\beta\in B}a_\beta(y)\chi_{u_\beta}(z).
$$
每个物理切片的非零二次卷积由整个零和组构成，不同旧差不可能合并。同一旧差的多个零和组可以合并，不影响总和。常数项为 $\sum_\beta a_\beta(y)^2=1$。因此每个切片的平方恒等于1，新函数严格布尔。每个真带只是标签改为 $u_\beta$，$S_\beta$ 全部保留，所以稀疏度恰相等，且原全谱线性秩不会下降。

任选 $\alpha_\beta\in S_\beta$，定义 $K=\operatorname{span}\{\alpha+\alpha':\alpha,\alpha'\in S_\beta,\ \beta\in B\}$，$\bar A:E\to Y^*/K$ 为 $e_\beta\mapsto\alpha_\beta+K$。其定义不依赖所选代表。新谱由所有 $(\alpha,u_\beta)$ 构成；向外部做同带差分得 $K\oplus0$。模去此空间之后，全部谱点张成为映射 $v\mapsto(\bar A v,v+R)$ 的像，核为 $R\cap\ker\bar A$，故真实全谱线性秩精确为
$$
 \dim K+|B|-\dim(R\cap\ker\bar A).
$$
允许所有零和分组时，最大全谱秩由上式中交维数的最小值给出。确实，任一合法联合标签提升 $e_\beta\mapsto u'_\beta$（要求存在投影 $u'_\beta\mapsto\beta$）按每个切片中实际的提升差分组，都得到零和划分，且 $R\subseteq\ker(e_\beta\mapsto u'_\beta)$。此 $R$ 的通用提升能够投回所给联合提升，秩至少一样大。反向每个构造已实际合法。

此公式是第5轮条件在指定商层上的版本；额外的外部自由度由 $\ker\bar A$ 精确处理。若有一般 A28 的零比值反例，可在任何指定的带分解上这样提升，保持 $s$ 并提高或保持谱秩；结论不依赖把某个非仿射活动事件当作物理仿射限制。它本身没有给出所需的普适交维数下界。

例如在四维 $B=\mathbb F_2^4$ 上，让两种频域符号相位是 $q_0(\beta)=\beta_0\beta_1+\beta_2\beta_3$，$q_1=q_0+\beta_0\beta_2$，系数均为对应符号乘 $1/4$。两者极化非退化，逆变换皆为布尔相位。用一个外部物理位 $y$ 作选择，则每个外带单点 $\alpha_\beta=\beta_0\beta_2$，$K=0$。取各非零差的所有边为一组，$R$ 为完整四点仿射关系空间，维11：模关系后所有基标签由0和四坐标点生成，而 $e_\beta\mapsto(1,\beta)$ 已给像维5。$\bar A$ 在 $R$ 上非零，因 $\beta_0\beta_2$ 在 $0,e_0,e_2,e_0+e_2$ 有奇校验。这样交维数是10，得出的全谱秩为6，直接相加商标签秩及同带差秩只给5。此例核验的是指定分组的公式，而非该输入所有分组的最大值。
