# 一条六边形边放宽后的连续输入宿主

在课题根目录执行下列命令，输出写入已有的可写目录 `$WORK`。基准凸体、曲线和法向取自同一版本的 computation。复算使用 Python 3.13.18 和 python-flint 0.9.0，须启用断言。两个外包回执的文件名用于随后的真实面积比较，归档同名副本放在 H 中。

```sh
H=computation/released-hexagon
python -I -B "$H/certify_released_chain_host.py" computation "$H/continuous-153-input.json" "$WORK/continuous-153-area.json"
python -I -B "$H/certify_released_chain_host.py" computation "$H/single-H1-input.json" "$WORK/single-H1-area.json"
cmp "$WORK/continuous-153-area.json" "$H/continuous-153-area.json"
cmp "$WORK/single-H1-area.json" "$H/single-H1-area.json"
python -I -B "$H/compare.py" computation computation/single-repair "$H/continuous-153-input.json" "$H/single-H1-input.json" "$WORK/difference.json"
python -I -B "$H/certify_counterexample.py" computation "$H/single-H1-input.json" "$H/counterexample.json" "$WORK/missing-box.json"
python -I -B "$H/certify_monotone_cost.py" computation "$H/single-H1-input.json" "$H/cost-input.json" "$WORK/cost-86.json" --gain 8.6e-8 --steps 320
```

每个构造均包含原来的连续凸核及150个已放好的定宽输入。对 \(w=1.2\) 的闭参数区间，两个端点的上包加上中间姿态二阶插值误差；单面方案按第0、2面固定平移。比较程序用4,096个弧段重建真弦和旧边，单面宿主、交错三面宿主分别使用22,408与22,405个法向。单面宿主的一个真外三角形给出额外下端。对每一个候选重算有向支撑和整数凸包，用一方的下端减另一方的上端给出严格差值。旧六面的端点强制链见第219轮。

遗漏闭盒的证书在每个闭角分支中按六向约束界住真实平移，再检验五个指定分离方向。增量证书对两个手性共计893个角叶分支，结合 AH 真弦外的实际材料与已知凸帽面积，给出相对同一宿主的增量。宿主的精确面积下端来自 difference.json。允许平面中任意方向扩张时，另一个定宽输入及整个参数矩形的上下界见[自由扩张证明](../../proofs/free-face-expansion.md)。
