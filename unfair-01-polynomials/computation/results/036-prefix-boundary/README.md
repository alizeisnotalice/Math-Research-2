# q21剩余族的两个独立前缀

正向 s8,t20：0<a≤1/2；倒数 s12,t20：0<a<1。完整量词与失败尾部见 proofs/036-paired-three-point-and-boundary.md。
两前缀没有拼接成同一有限因子。a=1/2时正向第44项5/8、反向第48项5/16，反向第44项已为1而非零尾H，故不是反例。

```sh
/Users/zzh/.venvs/main/bin/python verify_prefix_intervals.py --forward k2-prefix.json --reverse k2-reverse-prefix.json --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/prefix-replay036.json
```

该程序从卷积重建参数多项式，以有理Bernstein系数覆盖整个区间；本轮0.01895秒。verify_k2_prefixes.py另以Fraction重算a=1/2，回执已归档；其输出在脚本目录，复跑请先复制本包到CodexTemp。audit_prefix.wl在Mathematica15.0.0实际Reduce，并显式Exit[]；结果见 wolfram-prefix-audit.json。中间误把“最早零尾差异”与“最早非二元项”混同已在验证器修正，不沿用旧失败结果。
