# AH 连续路径、最短性与参数图像独立审计

审计基线：`7ae44feb6a38e9e69207b234d494af307138f546`。
本文件只审计第18--20轮固定 \(\sigma=1.5494^\circ\) 的 AH 局部构造；不把它解释为 Gibbs 全局万能性、面积认证或 Lebesgue 常数的新界。

## 严格结论：通过

1. 第18轮以有理区间构造的 \(Y,T,Z,C,G,E_3,N_-\) 及路径提升可重放。圆弧和切线段的四圆不等式、中心根的单值性与连续性均由严格余量承担。
2. 第19轮的外圆障碍等式
   \[
   |X-C|^2-1=\ell\{2(X-L)\cdot e_L+\ell\}
   \]
   与其37个闭参数盒的严格正性，足以推出整个参数图像位于 \(C\) 的闭单位圆外。证书的闭盒面积和为1，盒内且无内部重叠；由于有限闭盒并集在单位方形中闭，面积相等排除非空相对开补集，故确实覆盖整个方形。
3. 外圆最短性的势函数论证是自足的。对
   \[
   g(r)=\sqrt{r^2-1}-\arccos(1/r),\qquad F_\pm=\pm\phi+g(r),
   \]
   有 \(g'(r)=\sqrt{r^2-1}/r\) 及 \(|\nabla F_\pm|=1\)；按极角提升的绕数分两类，给出任何圆外曲线的长度下界。它不要求 \(R\) 预先是 Jordan 域。
4. 第20轮的切点不等式在唯一的 \((L,N)=(G,N_-)\) 角允许等号，且该等号另由解析恒等式处理。对于固定 \(X\)，\((X-L)\cdot e_L\ge0\) 和 \((X-N)\cdot e_N>0\) 分别选出各中心线上的唯一根；两根合并时仍唯一。因此连续参数映射是单射。紧矩形到平面 Hausdorff 空间的连续单射是到像的同胚，其边界像为 Jordan 曲线。
5. 原路径和配对路径的开参数象限分离仅在公共端点 \(Z\) 相交。文中所说的固定 \(s=\pm\theta\) 圆弧重合是路径遍历的往返段，删除后 Green 边界积分抵消的判断合理。

## 实际重放

在本独立共享对象克隆中，以下脚本均以 Python 3.15.0 实际退出0：

```text
computation/certify_ah_path.py
computation/certify_ah_shortest.py
computation/verify_ah_shortest.py
computation/certify_ah_coordinates.py
```

重放输出在 `/Volumes/data/CodexTemp/T2/lebesgue-universal-cover/bridge-next/`。最短性证书有37个叶盒，独立验证器复核37个盒、总面积1及所有两两内部不交。

`check_ah_high_precision.py` 未运行：当前解释器没有安装 `mpmath`。该脚本的80位输出在原设计中是数值诊断，不用于上述严格结论。

## 范围与仍需连接的部分

- 本审计没有独立验证 AH 片段与 Gibbs Definition 8 全部外部边界标记的逐点对应；现有第20轮证明的是当前坐标模型中该片段和其配对的拓扑。
- 本审计没有认证 Green 积分的严格面积误差、其余 EH 切除、去除后凸性或 orbiform 反射放置。
- 所有“重合段消去”只消去同一路径的反向遍历，不能替代对未审计边界段的简单性证明。

## 固定文件哈希

| 文件 | SHA-256 |
| --- | --- |
| `proofs/018-ah-continuous-path.md` | `54b0dc2fefa6459cf103bf603ce3f45a36e9444ec47df011b9b7daac736d5400` |
| `proofs/019-ah-shortest.md` | `6e764c3b7d9a40d7cc0c7755a260f9e728174cfd8cf11c2acde5aca720be99ca` |
| `proofs/020-ah-coordinate-topology.md` | `6a1621d8ff48cb05bba13a4931b2fa50efbc7210d1ce759efd36d2d08bf453c3` |
| `computation/certify_ah_path.py` | `98f50387de75d0a331f053fbb8d194f738fd8c68cfc0dcfec2915a94eb192b9a` |
| `computation/certify_ah_shortest.py` | `5211461f1a3fecf1f579802ec2c08080a29d699e05db6ef1100a93bd5c1e8c4e` |
| `computation/verify_ah_shortest.py` | `1f35ba1bbe663804204a84ce50c630c34bdfaa558e4cf99ffca68303a41589e7` |
| `computation/certify_ah_coordinates.py` | `5d89d49dfb03309c6c2b7e992613c84f11dee1e61dddbf86097b8b2a10fe8ef5` |
