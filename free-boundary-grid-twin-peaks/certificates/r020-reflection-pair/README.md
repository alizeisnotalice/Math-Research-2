# 第20轮：固定镜像标签对的有限反例

`reflection-pair-orbit-audit.md`给出4×4一个固定无序标签对轨道，全部256种朝向中恰双峰同侧12个、异侧6个；这只否定逐轨道same→split单射。`spin-encoding-audit.md`给出局部spin判据并以直接峰判定交叉核对。`near_central_split_witness.py`生成中央异侧而距离仅2的标号，提醒中央split事件不自动意味着线性距离。

三个Python脚本直接向stdout输出JSON；主Agent实际重跑spin判据及near-split构造，解析JSON与归档完全相等。执行脚本时将stdout重定向到 `/Volumes/data/CodexTemp/T2` 的独立子目录。没有原题反例或渐近推论。
