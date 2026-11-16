# R8–L9 的有理根圆盘证书

对象为 R=1+x+x²+x³+x⁴+x⁵+x⁸ 及 C=R(x)R(x⁹)。证书仅认证 C 的80个互异非实根：8个母根与72个提升根，按共轭组成40组。它不认证全部根组分配的系数符号。

十进制根建议转成精确有理中心 c；所有半径 r 严格正。母根使用
|R(c)|+Σ(k≥2)|R^(k)(c)/k!| r^k < |R'(c)|r。
提升根使用 z⁹−α，其中 α 位于已认证母圆盘；左边额外加母圆盘半径。各不等式由 Rouché 定理给恰一个根。逐对平方距离严格大于半径和平方，结合次数8及每母根9次提升，保证完整覆盖。

模长界使用 Fraction 和整数 isqrt：若 k=floor(s|z|)，则 k/s≤|z|<(k+1)/s。代码以 k=isqrt(floor(s²|z|²)) 精确实现，s=10^160。验证器重算全部界、目标整数卷积、输入绑定、覆盖及分离；不信任证书success或slack。生成器与验证器共享 root_ball_common.py，主线程另行审查数学论证，不能称为完全独立双实现。

从仓库根目录复现（输出在data临时目录）：

```sh
PY=/Users/zzh/.venvs/main/bin/python
BASE=unfair-01-polynomials
WORK=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round022-parent
PYTHONDONTWRITEBYTECODE=1 $PY $BASE/computation/root_balls/verify_root_balls.py --certificate $BASE/computation/results/022-R8-L9-root-balls-certificate.json --input $BASE/computation/results/021-R8-L9-inputs.json --record $WORK/recheck-root-balls.json
PYTHONDONTWRITEBYTECODE=1 $PY $BASE/computation/root_balls/test_negative_cases.py --certificate $BASE/computation/results/022-R8-L9-root-balls-certificate.json --input $BASE/computation/results/021-R8-L9-inputs.json --record $WORK/recheck-negative.json --work $WORK
```

负例实际覆盖负/零半径、重写所有绑定哈希后的错误目标、改变中心、重复来源覆盖，全部被拒。未安装或使用 Arb/python-flint；其在本轮环境中不可导入。
