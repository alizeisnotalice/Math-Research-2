# 前五行局部排除的代数整数反例

令ω=(1+√5)/2，`input.json`中每项整数对(m,n)表示m+nω。两侧非零系数均为代数整数且严格在(0,1)，常数项1，A5=B5=0，乘积C0,…,C5全部1。但是

C6=−31770+14208√5∈(0,1)，

因此它只否定“仅前五行二元预算加A5=B5=0即矛盾”的过强局部引理。它不是原题反例，也不说明全模型的(r,t)=(1,5)可行。

`generator`保留SymPy1.14.0实际生成与核验的v2输出；旧√3/2例不满足代数整性，原运行未改写，v2另起。独立`verify_fraction.py`不导入生成器，用完整O_K整数乘法重新计算全部卷积、首一整最小多项式及√5隔离区间平方证。其实际Fraction回执为`fraction-audit.json`。

```sh
python verify_fraction.py --input input.json --out /Volumes/data/CodexTemp/T2/new-local-check.json
```

所有严格正性和小于1界均为有理数比较；没有浮点判断。第六项非二元是第一处失败，完整乘积不能冒作均匀和分布。
