# 完整周期水平集的准确计算

仅需 Python 标准库（Python ≥3.9）。所有原始坐标和质量使用有理字符串；浮点输出只用于显示。`periodic_halo.py` 计算严格内界和完整外界；`independent_periodic_halo.py` 独立计算全部盒质量及矩形并面积。

从本课题目录运行，输出放 data 临时目录：

```sh
python3 -I -B computation/periodic_halo.py computation/inputs/legacy_r004_transfer_34_27.json --bits 40 --tail-radius 5 --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/primary.json --rectangles-out /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/rectangles.json
python3 -I -B computation/independent_periodic_halo.py --input computation/inputs/legacy_r004_transfer_34_27.json --bits 40 --radius-cap 5 --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/independent.json --compare-primary /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/primary.json --compare-rectangles /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/rectangles.json
python3 -I -B computation/tests/test_periodic_halo.py
python3 -I -B computation/verify_saved_results.py
wolframscript -file computation/symbolic_checks.wl
```

运行前确认输出目录已创建、data 已挂载，并把 TMPDIR/TMP/TEMP 及缓存指向 data。Mathematica 脚本结束显式退出。

前两条命令重新枚举几何；`verify_saved_results.py` 只核对归档结果之间的准确算术和文件身份，不重算几何。原始矩形大数组可重建，不进入公开提交。软件版本和实际参数分别记录在结果 JSON 中；不同算法使用不同 Python 版本的事实按实际保留。

一般 L¹ 下界还依赖 [完整解析传递](../proofs/periodic-complete-halo.md)，不能只凭程序输出宣称。结果中的周期比值外界不适用于全部输入，也不是 C₂ 的一般上界。
