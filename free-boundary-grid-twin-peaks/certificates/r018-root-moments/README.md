# 根矩的整数证书

`root_moments_ie.cpp` 复用既有独立集计数核，在 n≤6 完整独立集上累计新矩；`n7_root_moments_d4.cpp` 使用 D4 轨道进行另一累计。它们的计数递推来源相同，不声称是彼此独立的数学方法。n4..6 全表与轨道的四个整数一致；n7 A2/D 与前轮距离证书一致。

```sh
clang++ -O3 -std=c++17 -pthread root_moments_ie.cpp -o /Volumes/data/CodexTemp/T2/a022-r018-moments
clang++ -O3 -std=c++17 -pthread n7_root_moments_d4.cpp -o /Volumes/data/CodexTemp/T2/a022-r018-d4
/Volumes/data/CodexTemp/T2/a022-r018-moments 6
/Volumes/data/CodexTemp/T2/a022-r018-d4 7
python3 verify_root_moments.py
python3 verify_d4_root_moments.py
```

两个 Python 脚本在本目录核验既有 JSON 并更新小摘要。C++ 输出新 JSON 到 stdout；资源文件来自实际已完成运行。n7 最大 RSS 为 2,833,793,024 bytes；4,500,017,160 bytes 是峰 memory footprint，两者不得混用。主审又编译归档版并复算 n4，见 archive-root-checks.json。证明与无穷子列条件见 `../../proofs/root-moment-inclusion-exclusion.md`。

本轮只研究自由边界方格、原图距离、均匀标号条件 K=2。有限数据不证明渐近结论。

复现环境为 Python 3.15、Apple Clang 21、Swift 6.4.2 / Metal。所有编译与输出请放 `/Volumes/data/CodexTemp/T2`，Clang/Swift 模块缓存使用 `/Volumes/data/CodexTemp/cache`。
