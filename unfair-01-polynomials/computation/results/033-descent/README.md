# 短因子正性下降复核包

解析结论见 [完整证明](../../../proofs/033-short-factor-descent.md)。交叉替换假设下，当互素 \(p<q<2p\) 且短母因子次数小于 \(\min(2q,3p)\) 时，两母因子均非负；另有 \(p=3,q=5\)、短侧次数小于10的更强结论。另一侧次数无上限。一般P31与均匀和猜想仍开放。

independent-final-review.md 绑定最终证明SHA，审查缩放前后的首一性、常数情形、倒数反转、首负/末负配对及完整边界。binding.json 绑定正文与实际输出。

主线程在现有 Python 3.13.18 数学环境执行 check_short_factor_identities.py，以整数稠密卷积检查54,528个系数公式，并检查各被覆盖参数的碰撞指数。最终量词强化后重新运行一次，耗时约0.38秒，进程退出0。有限索引/系数核对只是检查实现和公式，普遍结论由解析证明承担。

    /Users/zzh/.venvs/main/bin/python check_short_factor_identities.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/short-factor-recheck.json
