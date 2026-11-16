# Hd(3,4) 全二元参数宏行证明的精确核验

证明见 ../../../proofs/027-Hd-x3-x4-all-binary.md。d≥6目标二元恰排除6/8/9/12。d7由本轮完整根证书，d10由第25轮解析证明；其余全部d≥11、d≠12由本轮统一宏行证明。

在本目录运行：
```sh
python3 check_hd34_macro_carries.py
python3 Hd34_tail_independent_audit.py
```

实际Python3.13.18、SymPy1.14.0；主线程两脚本复跑均退出0，回执parent-audit-runs.json。第一个脚本会重写同目录的可再生有限账本JSON；第二个从头重建支持，不读第一个脚本/结果。它们核验完整有限carry、W消失边界、行恒等式与周期36的gcd；一般d的推导与UFD由证明正文承担，不将有限扫描本身当一般证明。
