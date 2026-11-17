# r030 稳定纤维长度—峰距联合计数

核心定义、严格有限结果、解析范围、资源和复现命令见 `REPORT.md`。

- `gap_geometry_dp.cpp`：精确 \((L,d,e^2)\) 两阶段 subset DP。
- `gap-n3-all.json`、`audit_n3_joint.py`、`n3-joint-audit.json`：独立全排列核验。
- `gap-n5-v*.json`：六个 D4 孔类结果。
- `summarize_n5_geometry.py`、`n5-gap-geometry-summary.json`：D4 合并与 r018 距离/矩交叉核验。

本目录是临时研究产物；未修改共享 Git，也不包含任何渐近外推。
