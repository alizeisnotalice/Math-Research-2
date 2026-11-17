# 第 4 轮：有限 Walsh 符号码的子空间审计

`signcode_subspaces.py` 穷尽所有 \(n\le4\) 的布尔函数，按完整的坐标绝对 Walsh 向量分组。每一组以固定代表元建立相对符号集合 \(K\)，并从整体符号及输入平移建立子码 \(C\)。程序精确检查 \(C\subseteq K\)、\(K+C=K\)，再在 \(K/C\) 中穷尽所有全体元素仍属于该集合的线性子空间。

对每个最大见证，记录频率支撑、未归一化整数 Walsh 系数、\(C\) 基与补空间符号码基。独立验证器 `verify_signcode_witnesses.py` 对 2299 个组的全部见证逐个构造图提升的每个切片并做整数逆 Walsh；2302 个切片均为 \(\{-1,1\}\) 真值表。

MacBook 的实际运行目录是
`/Users/zhengzhihao/CodexRemoteMath/a28-selector/r004`，命令为：

```sh
/usr/bin/time -l /usr/bin/python3 signcode_subspaces.py \
  --max-n 4 --out remote-result.json
```

该完整扫描用时 `1.02` 秒（用户 CPU `0.91` 秒），最大 RSS `55721984` 字节。源 SHA-256：
`0b5079dfd26a90ef1a6ac41d781e806b2a07e6609d7ee7aa39f67ef32566839e`；结果 SHA-256：
`50166ce3632c799494a889d5de88ad1e5eef409ee940409f31ee72ebf2cb6f47`。

扫描中最大商维数为 11；相应组的商集合仅含 28 元素，DFS 实际检查 73 个子空间。按支撑大小，最大见证为 \(s=16\)、\(n=4\)、\(\dim C=5\)、额外维数 \(r=2\)、\(\dim V=7\)。因此有限范围内，符号切换确实可超出平移和整体符号产生的秩，但这一现象本身不证明任意维结论。

远端结果在
`/Volumes/data/CodexTemp/T2/a28/r004/remote-result.json`，本机独立收据在
`/Volumes/data/CodexTemp/T2/a28/r004/local-witness-receipt.json`。运行结束后远端没有本任务进程。

主Agent归档补充：MacBook使用Python 3.9.6。本目录保存取回的原始结果及上述收据。这里的dimV是谱点齐次提升秩；n=4、s=16见证的默认图提升只有6个变量，实际线性秩和仿射维数均为6，并非线性秩7。另乘一个独立字符后线性秩才为7。主Agent的`verify_bent_quotient.py`实际核验了这个区别，独立逆谱穷举全部65536个符号并排除三维商子空间，结果见`parent-bent-verification.json`。

`triangular_lifts.py`构造任意k的三角置换提升；本轮实际运行k=0..5，逐个非零差值检查所有无序谱点对乘积和为0，最大s=1024、d=63。k≤3另作完整逆Walsh；置换参数k≤3穷举，k=4,5各取固定种子20261116+k的2000个样本。一般任意k结论来自正文证明，样本不是一般证明。结果见`triangular-results.json`。两项本机程序使用Python 3.15.0，整数计算、零容差，均正常退出。
