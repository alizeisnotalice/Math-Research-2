# 第 3 轮：新选择器节点树的有限精确核验

脚本 `verify_fresh_selector_trees.py` 对叶标号
\(\pm\chi_b\) 的新选择器二叉树直接构造 \(\{-1,1\}\) 真值表，使用整数 Walsh–Hadamard 变换取得 Fourier 支撑，并对支撑向量的二元线性差与齐次提升向量做消元。没有浮点阈值。

检验内容包括每个内部节点的：

- \(4s\ge(\operatorname{affdim}S+2)^2\)；
- 以 \(A_i=\operatorname{span}\{(\alpha,1):\alpha\in S_i\}\) 给出的非例外递推
  \(a_H=a_L+a_R-h+1\)；
- 纯数据共同频率数 \(M\) 的精确损失界
  \(2M\le2s_L+2s_R-s_H\le3M\)；
- 删去递推常数 `+1` 的负面测试。

例外定义为两个子函数都不含其私有选择器频率；此时显式保留常量、相同字符、相反字符和根平移字符等边界例。

MacBook 上的实际完整运行目录为
`/Users/zhengzhihao/CodexRemoteMath/a28-selector/r003`，运行命令为：

```sh
/usr/bin/time -l /usr/bin/python3 verify_fresh_selector_trees.py \
  --depth3-samples 2000 --out remote-result.json
```

使用 Python 3.9.6，实耗 `16.03` 秒，用户 CPU `15.91` 秒，最大 RSS `17137664` 字节；本任务进程已退出。源 SHA-256：
`41a6fa7cef00bd9441c8448d336f4c03a8d54de0c5e22db2ea169122a3adeb92`。
结果 SHA-256：
`7cd57133c27810979e92bc562338252ae256bf119231e6a9a730957430944643`。

精确覆盖是四叶平衡树 4096 标号、四叶梳形树 4096 标号；八叶平衡树为种子 `20261116` 的 2000 个样本。结果位于
`/Volumes/data/CodexTemp/T2/a28/r003/remote-result.json`。这是有限模型核验，不构成一般定理的证明。
