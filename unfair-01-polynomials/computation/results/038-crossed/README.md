# 双侧整数切片下降与辅助命题精确反例

完整数学论证见 `proofs/038-integer-slices-and-signed-mothers.md`（相对课题根目录）。主线程的 `check_signed_cross.py` 只用稀疏整数卷积：核验79个互素参数的普遍带符号母式族，以及p2/q3的完整带进位例。独立 `check_two_sided_projection.py` 使用SymPy1.14核验相同结构及非互素p4/q6整数递推。

复现：在本目录运行 `python3 check_signed_cross.py` 与 `python3 check_two_sided_projection.py`。脚本结果写同名JSON；主线程已从本归档代码实际执行。`search_carry_counterexample.py` 保留找到例子的有限整数探索参数；它不宣称最小性。

P31“完整交叉非负迫使母因子非负”已被该族精确否定，包括不同奇素数p3/q5；原交叉因子A/B全部二元，故原均匀和猜想没有被否定。双侧零余类整数性仍是充分下降条件，不是已证的一般性质。
