# 前缀移动核的有限核验

```sh
python3 audit_resample_move.py
clang++ -O3 -std=c++17 prefix_connectivity.cpp -o /Volumes/data/CodexTemp/T2/a022-r018-prefix
clang++ -O3 -std=c++17 two_step_endpoint_connectivity.cpp -o /Volumes/data/CodexTemp/T2/a022-r018-two-step
/Volumes/data/CodexTemp/T2/a022-r018-prefix
/Volumes/data/CodexTemp/T2/a022-r018-two-step
```

n2/n3 全部84个 (n,r,tau) 支持实例完整枚举。两步端点核的中间态允许 birth 无效，只要求两个端点有效；这对保持 Q² 对称性是必要的。小图连通性不推出任意 n 的不可约性或混合速度。Python 粒子校准使用其自身浮点实现，只是诊断，不等于 C++/Metal 整数选父认证。

本轮只研究自由边界方格、原图距离、均匀标号条件 K=2。有限数据不证明渐近结论。

复现环境为 Python 3.15、Apple Clang 21、Swift 6.4.2 / Metal。所有编译与输出请放 `/Volumes/data/CodexTemp/T2`，Clang/Swift 模块缓存使用 `/Volumes/data/CodexTemp/cache`。
