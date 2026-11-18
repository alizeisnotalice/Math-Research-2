# 两个固定宿主的可放性和面积

H155 是第155轮圆盘、四个正定宽奇边形及三十点的明确共同凸包。H175 是第175轮七输入的明确共同凸包。两个来源姿态和解析支撑的源文件由 source-manifest.json 锁定。certificate-labels.json 给出原角格的三个有理单位法向编号；五点输入有少量重新选择的三元组。所有权重都在重放时从行列式精确恢复，并验两坐标的加权和为零。每个输入各有两种手性。

在可写的项目副本的课题目录运行，Python 3.13.18、python-flint 0.9.0；脚本使用断言。角格重放工作精度为160位，所导入的历史几何模块保留自身初始化产生的区间误差。

```sh
C=computation/two-host-witnesses
python -I -B "$C/smallest/verify_four_point_reduction.py" computation "$C/smallest/four-point-reduction.json"
python -I -B "$C/area_cert.py" computation
python -I -B "$C/point60/replay.py" computation --reuse-checked
python -I -B "$C/certify_secants.py" computation "$C/secants.json"
python -I -B "$C/derive_gains.py"
python -I -B computation/rotation-modulus/certify_rotation_moduli.py computation/rotation-modulus/rotation-moduli.json
```

第二条命令从几何输入重建旧角格估计和四项明确摆放的面积上端。第三条命令用已独立核过的八项全角回执并验证其中全部下端，重建两个35点集。省略 --reuse-checked 可同时再跑全部八项。第四条命令用同一组三方向的双端点公式从头核对16种宿主、输入与手性组合，先对五点与十一点做有理距离检查。末步面积转移采用内圆与支撑线外点相切的三角形；新增面积均指包含旧宿主的自由凸扩张。部分旧回执保存实际运行时间，复跑时该字段随机器负载变化。
