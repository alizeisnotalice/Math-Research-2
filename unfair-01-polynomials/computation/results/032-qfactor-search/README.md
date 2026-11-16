# P31的有界有理因子探索

实际范围p3q5、18≤degH≤24；H首一二元、H(0)=1且H(x³)H(x⁵)二元。完整支撑DFS得到19,511个H；按SymPy有理不可约因子及重数枚举68,996个分配。28,208个带符号分配每个都有精确整数负系数见证，没有双交叉非负候选。18–20次现已被本轮解析定理覆盖；新的探索是21–24次。

**范围限制：** 不可约因子清单的完备性仍依赖SymPy1.14.0的QQ分解后端。搜索和独立验证程序分离、逐项重新卷积，但二者使用同一分解后端；未把后端重复同意当作外部不可约性证书。这里是可复核的有理输入探索，不能作为任意实因子或一般P31的排除结论。独立验证严格确认的是所存分配与负系数见证。

逐H因子列表、重数、全部分配编号和负项保存于四个确定性gzip账本。`compression-manifest.json` 同时保存压缩与原始SHA、体积；压缩不改变证据。实际Python3.13.18、SymPy1.14.0。脚本不使用随机数；单CPU搜索与验证均已退出。

复现（从本目录；已有数学环境）：

```sh
python3 unpack_ledgers.py --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round032-reproduce
/Users/zzh/.venvs/main/bin/python verify_p31_qfactor.py --search-out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round032-reproduce/search-d20-d24 --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round032-reproduce/verification-20-24.json
/Users/zzh/.venvs/main/bin/python verify_p31_qfactor.py --search-out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round032-reproduce/search-d18-d19 --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round032-reproduce/verification-18-19.json
```

全量重新搜索的参数见 `search_p31_qfactor.py --help`，应选择新的临时输出目录。DFS每次添加t的完备性来自：若旧支撑已有两种和表示，添加点不会消除碰撞；反之所有无碰撞支撑的每个递增前缀都会保留。新点的2|S|+1个和须互异且不与旧和重叠。不能把超时、漏项或后端未知写成排除。
