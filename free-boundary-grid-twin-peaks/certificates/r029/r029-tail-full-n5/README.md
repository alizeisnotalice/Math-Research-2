# r029 穿孔首邻点完整尾分布

核心结论、范围、运行资源和复现命令见 `REPORT.md`。本目录只含临时研究材料，未改动共享 Git。

- `tail_full_dp.cpp`：精确 C++ 两阶段 subset DP。
- `tail-full-n3.json`、`tail-full-n4.json`：对 r028 完整直方图的独立交叉核验。
- `tail-full-n5-*.json`：六个 D4 孔类的精确 \(n=5\) 结果。
- `n5-tail-summary.json`：D4 合并后的完整全局尾。
- `audit_tail_full_dp.py`、`summarize_n5_d4.py`：有理数审计及合并程序。
