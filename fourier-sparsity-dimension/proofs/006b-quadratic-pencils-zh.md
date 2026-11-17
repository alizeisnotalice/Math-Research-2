# 第6轮补充：二次相位子类的范围和包含极大性

对2k变量的二次bent相位f₀与二次掩码线性空间V，设f₀+V全部bent且V含全部仿射函数。ANF环境维数直接给dimV≤1+2k+binom(2k,2)=2k²+k+1。因此该子类只有O((log s)²)维，不能产生超平方根谱秩。此粗界无需矩阵空间的外部定理。

极化将二次相位模仿射函数同构到2k阶交替矩阵；二次相位bent当且仅当其极化非退化。取下述块三角矩阵族，其方向维数k(k−1)，加全部仿射项得到dimV=k²+k+1。

下述一般k的证明只说明这一特定构造在二次子类中按包含关系极大，不能排除别的更大二次族；它也不排除添加三次或更高次掩码。特别是k=3时该构造维数13，小于第4轮三角MM的15。

# 标准二次 symplectic pencil 的包含极大性

令 \(X,Y\cong\mathbb F_2^k\)，并以块矩阵写交替形式。令
\[
 \mathcal P=\left\{
 \begin{pmatrix} A&I+N\\I+N^T&0\end{pmatrix}:
 A\in\operatorname{Alt}_k(\mathbb F_2),\ N\in\operatorname{NT}_k(\mathbb F_2)
 \right\},
\]
其中 \(\operatorname{NT}_k\) 是严格上三角矩阵空间。每个元素非退化，且方向
空间维数为 \(2\binom k2=k(k-1)\)。

**命题。** 若交替矩阵
\[
D=\begin{pmatrix}C&X\\X^T&Z\end{pmatrix}
\]
满足 \(P+D\) 对每个 \(P\in\mathcal P\) 都非退化，则
\(Z=0\)、\(X\in\operatorname{NT}_k\)，故 \(D\) 已在 \(\mathcal P\) 的方向空间中。
特别地，\(\mathcal P\) 没有保持全体非退化性的真一维线性扩张。

**证明。** 写 \(U=I+N\)、\(K=A+C\)。因 \(A\) 任意交替，\(K\) 也任意交替。

若对某个 \(U\)，\(P_0=U+X\) 奇异，取 \(K=0\)。则
\[
 \begin{pmatrix}0&P_0\\P_0^T&Z\end{pmatrix}
\]
被 \((u,0)\) 湮灭，其中 \(0\ne u\in\ker P_0^T\)，矛盾。因此在假设下每个
\(U+X\) 必可逆。

若 \(Z\ne0\)，取 \(v\) 使 \(Zv\ne0\)，固定任一可逆 \(P_0=U+X\)，并令
\(u=P_0^{-T}Zv\)。则 \(u\ne0\) 且
\[
 u^TP_0v=v^TZ^Tv=v^TZv=0,
\]
因为 \(Z\) 交替。故存在交替 \(K\) 使 \(Ku=P_0v\)：映射
\(K\mapsto Ku\) 从 \(\operatorname{Alt}_k\) 满射到 \(u^\perp\)。于是
\[
 \begin{pmatrix}K&P_0\\P_0^T&Z\end{pmatrix}
 \binom{u}{v}=0,
\]
又矛盾。所以 \(Z=0\)。

于是 \(Y:=I+X\) 满足 \(Y+N\in\mathrm{GL}_k(\mathbb F_2)\) 对所有严格上三角
\(N\) 成立。特别地 \(Y\) 可逆，并对任意 \(i<j\)，\(Y+E_{ij}\) 可逆。矩阵行列式引理在
\(\mathbb F_2\) 给出
\[
 1=\det(Y+E_{ij})/\det Y=1+(Y^{-1})_{ji},
\]
所以 \(Y^{-1}\) 的全部严格下三角项为零。故 \(Y\) 为上三角；其每个对角元非零，只能为1，
即 \(Y=I+N_0\) 且 \(N_0\) 严格上三角。因此 \(X=N_0\)。\(C\) 本来任意交替，正是原方向空间
的左上块，完成证明。

范围：这是一个特定 pencil 的包含极大性，不排除另一、维数更高的 \(\mathbb F_2\) 仿射
symplectic 空间，因而不能替代小域中的全局 sharp 上界。


## 计算核对

本轮Python 3.15.0标准库实际枚举Alt₄(F₂)全部64个矩阵，其中28个非退化；所有仿射子空间按方向维数0、1、2分别有28、378、315个全非退化者，维数3以上为零。因而k=2时最大二次相位族维数为2+5=7。

对k=3标准构造，全部64个矩阵非退化；枚举32768个交替方向并排除原方向空间64个，剩余32704个均不能加入。该有限检查支持上面的一般包含极大性证明，不是k=3的所有矩阵族分类。

脚本与实际输出见computation/round006/quadratic。原运行及源码哈希登记execution-record.json；复现命令为python3 alternating_affine_search.py --out <data临时JSON>及python3 standard_pencil_extension.py --out <data临时JSON>。

本轮检索到的更锐交替矩阵维数结论带有域大小假设，不能直接用于F₂。这里完全不用该锐结论，也不宣称F₂上k(k−1)是一般最大方向维数。
