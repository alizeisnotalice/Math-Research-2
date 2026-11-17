# Sprague 固定角切除审计

审计对象为 `audit-sprague/proof.md` 与对应有理区间程序，角度固定为
\(\sigma=1.5494^\circ\)。本审计不改来源文件，也不把固定角切除结果写成新的
Lebesgue 常数界。

## 结论

通过。未发现支持接触、fan 内部排除或闭包补点中的实质漏洞。

### 接触段

`P` 仍位于原正六边形中。对任一定宽一的 orbiform，原六边形任一对平行支持线
之间距离恰为一，支持条带论证强制在两条原面都有接触。C/E 截角分别给
\(h(m_5)\le1/2\)、\(h(m_3)\le1/2\)；因为 \(m_2=-m_5\)、\(m_0=-m_3\)，
定宽关系迫出对侧旋转截线的支撑值至少 \(1/2\)。原面接触严格处在该截线的另一侧，
凸性遂给出 \(F_3F_2\)、\(B_3B_2\) 的接触。这避免了把这两段误当作 \(P\) 的原面。

对 DC，审计正确补上了不能从 \([D_1,C_2]\) 接触直接推出 \([D_1,G]\) 接触的缺口：
若接触落在 \((G,C_2]\)，则其相对 \(F_3F_2\) 接触距离严格大于一，矛盾。

### 曲边 fan

每一条圆弧都由端点内积正性选为短弧。对外部 apex 与弧点 \(c+u\)，端点条件
\((A-c)\cdot u\ge1\) 通过短弧的连通球面帽推广到整弧；扇形内点的平方距离恒等式
给出严格在相应单位圆盘外。Voronoi 线性泛函对 radial-chord 参数的256个完整区间
给出严格负上界；接触端点的零值另用向内导数处理，没有将数值零当成严格不等式。

ES 与 AS 的两扇由 apex--圆交点对角线分开。三条严格定向叉积将对角线置于两个外
射线之间、总角小于 \(\pi\)，故两扇内部分离而只在该对角线相会。fan 的有理盒也将
CS、ES、AS 三个闭区域分离。因此这些是完整曲边区域，而非仅为图上有限顶点的检查。

### 闭包补点

令 \(U=C_S\cup E_S\cup A_S\)。每个闭 fan 是其内部的闭包，且三者分离，因此
\(U=\overline{\operatorname{int}U}\)。Propositions 3--5 先排除
\(K\cap\operatorname{int}U\)。若 \(\operatorname{int}K\) 遇到 \(U\)，其开放性会迫使其遇到
\(\operatorname{int}U\)，矛盾。对 \(z\in\operatorname{int}K\) 和任意 \(q\in K\)，
\((1-t)q+tz\in\operatorname{int}K\)（\(t>0\)）并趋于 \(q\)，所以 \(K\subset
\overline{P\setminus U}\)。这正确避免了错误恒等式
\(\overline{P\setminus U}=P\setminus\operatorname{int}U\)。

## 实际运行

在 `lebesgue-sprague` 工作区实际运行：

```sh
/Users/zzh/.venvs/main/bin/python -B audit-sprague/certify_sprague_fixed.py
```

退出状态为0。重放 JSON 记录五条弧、CS/ES/AS 三个 fan 盒，以及固定角的所有严格
半平面、距离、Voronoi 与根分支检查。

## 与第26轮底形证明的衔接

主克隆 [第26轮 Pál 基覆盖证明](/Volumes/data/CodexProjects/T2/.local/workspaces/lebesgue-main/lebesgue-universal-cover/proofs/026-base-universal-cover.md)
独立证明任意直径至多一集合可放进未进一步切除的 \(P(\sigma)\)。本审计不替代其定宽补全。
反过来，本审计仅在 Gibbs 原文 Proposition 2/Corollary 1 的 orbiform/反射归约范围内说明
前三个 Sprague fan 可去除；AH/EH 切除、全局反射放置、完整凸性和最终面积仍不由本文件认证。

## 固定哈希

| 文件 | SHA-256 |
| --- | --- |
| `audit-sprague/proof.md` | `528d473d8c94999da1cd3b61064588b4f3d0b49173b09f65073b5405cc2c6adf` |
| `audit-sprague/certify_sprague_fixed.py` | `ef4f588c63b5fbaebff59e24678274940939cb89a4e38a9d26105928c75023ed` |
| `audit-sprague/rational_trig.py` | `97e0076b911f51a9c532e14cacb9fb6e3ae7e68707a3c5af77048f9c6623e87f` |
| 主克隆 `proofs/026-base-universal-cover.md` | `51ecba172c597557fec2918c1d07d9157957cb8ebaf221fcce57da1ea5b3ad6a` |
