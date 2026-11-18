# R85 高横向秩割的精确检查

[解析证明](../../proofs/085-high-transverse.md)给出三个现象。贪心预算之比可趋于 $2/3$。真实十点种子在横向余维五与六可以取等，并且满仿射包含实际谱中缺失的高秩标签。完全二叉树进一步给任意高横向余维的实际取等族，所有实际带秩二而缺失标签可达秩六。

对每个真实布尔种子精确遍历旧齐次子空间，求出原始 $L$。随后在新地址横向余维 5–11 上比较最小贪心谱基预算及按真实旧模长加权的逐带预算。若两项不同则分别记录。实际父函数由真实地址相位及数据平移重建，独立整数 Walsh 变换验证相位支持和绝对系数。

| 运行 | 精确结果 |
| --- | --- |
| Mac mini CPU，12 进程 | 40,606 个投影割，5,726 个贪心预算不足，真实模长加权预算无不足 |
| MacBook CPU，15 进程 | 81,281 个投影割，11,564 个贪心预算不足，加权预算无不足，一次取等 |
| 固定十点种子、固定首个秩二块的完整子族 | 横向余维 5、6 的 59,648,512 个矩阵族与旧子空间配对。按实际各带秩分类压缩后逐类剔除加权预算不足 |
| 独立保秩扩展搜索 | 双机 91,190 个真谱矩阵族，其中 9,986 次求出兼容线性空间维数。5,143 个存在非局部保秩扩展，均无一轴反例 |
| 两台 Metal | 每台独立重算 18,003 个割、408,871 条实际带的秩和预算 |
| 两台 Metal 真实父 Walsh | 每台逐系数核对五个实际布尔父函数的 35,328 个有符号整数系数，两台摘要一致 |

独立保秩扩展确实构造了一个横向余维六的取等割。对这份父谱的全部 23,936 个横向余维至少五的投影组合逐一核对，低余维由第 84 轮定理处理。另一个显式横向余维六的例子不满足满仿射包秩二假设。完全二叉树的三个实例分别达到横向余维 8、16、32，实有谱点数 64、256、1,024。第一个完整重算真值表，后两个核对全部稀疏卷积并对各 512 个真实平移布尔值做独立核验。固定矩阵子族的阴性结论是穷举证据，随机阴性只覆盖指定样本。

MacBook 主扫描平均使用 10.29 个 CPU 核，Mac mini 为 5.94 个。Metal 割预算的真实命令区间为 M6 10.258 毫秒、M5 Pro 3.017 毫秒。完整有符号谱的 Metal 命令区间另为 11.342 与 14.847 毫秒。均由 `MTLCommandBuffer` 读取。两台机器的相关运算已退出。

源文件使用 `../round084/phase_lifts.py`。Metal 后端沿用 `../round077/check_odd_walsh_metal.py`，SHA-256 为 `5ef8e0baf5a142512614cbaea22511927f8262ce0cf7155e1095a363ee8a4f9a`。MacBook 以已认证 SSH 标准输入接收源码，结果写回 Mac mini。等号并列时所留的随机见证可能随进程完成顺序变化。对相同的逐谱哈希，使用已归档的 `signed-witnesses-input.json`。

```sh
PYTHONPATH=../round084 python3 -B transverse.py scan --workers 12 --draws 768 --domain 0 > transverse-mini.json
set -o pipefail
PYTHONPATH=../round084 python3 -B transverse.py bundle --phase-source ../round084/phase_lifts.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - scan --workers 15 --draws 1536 --domain 1 > transverse-macbook.json
PYTHONPATH=../round084 python3 -B verify_witnesses.py transverse-macbook.json > witnesses-mini.json
PYTHONPATH=../round084 python3 -B compatible.py --workers 12 --draws 640 ../round084/phase_lifts.py transverse.py transverse-macbook.json witnesses-mini.json > compatible-mini.json
PYTHONPATH=../round084 python3 -B compatible.py --bundle ../round084/phase_lifts.py transverse.py transverse-macbook.json witnesses-mini.json | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - --workers 15 --draws 1024 --domain 1 > compatible-macbook.json
PYTHONPATH=../round084 python3 -B audit_extension.py > extension-mini.json
PYTHONPATH=../round084 python3 -B tree_high.py > trees-high-mini.json
python3 -B transverse_metal.py --witness-file transverse-macbook.json ../round084/phase_lifts.py transverse.py ../round077/check_odd_walsh_metal.py > transverse-metal-mini.json
python3 -B transverse_metal.py --bundle --witness-file transverse-macbook.json ../round084/phase_lifts.py transverse.py ../round077/check_odd_walsh_metal.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > transverse-metal-macbook.json
python3 -B witness_metal.py ../round084/phase_lifts.py transverse.py verify_witnesses.py ../round077/check_odd_walsh_metal.py signed-witnesses-input.json > witness-metal-mini.json
python3 -B witness_metal.py --bundle ../round084/phase_lifts.py transverse.py verify_witnesses.py ../round077/check_odd_walsh_metal.py signed-witnesses-input.json | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > witness-metal-macbook.json
```

## 回路、子空间和精确证书

下列单机程序用 Python 3.15.0 从本目录运行，以 python3 -B 加程序名即可重建相应记录。六、七维矩阵型证书用整数重新生成并验证；核验稀疏父函数或有限矩阵型的结论均按各自输入范围陈述。

| 程序 | 输出 |
|---|---|
| check_components.py | component-results.json |
| check_rank_four.py | four-results.json |
| exact_rank6.py | rank6_certificates.json |
| exact_rank7.py | rank7_certificates.json |
| verify_schur.py | schur_verification.json |
| verify_records.py | records.json |
| verify_synthesis.py | synthesis.txt |
| check_identification.py | identification.json |
| check_petersen.py | petersen.json |

双线性最大实现核验使用两个不交分段，覆盖四标签、指定三维 Z 的全部 16⁴ 组行子空间。共 61,377 组张满，其中 35,389 组逐标签评价满射；逐例核对可赋值空间维数。另核对 k=2,4,6,8 的显式张满回路。

~~~sh
python3 -B circuit_fields.py 12 0 7
ssh -o BatchMode=yes macbook-math /opt/homebrew/bin/python3.15 - 15 7 16 < circuit_fields.py
~~~

对应输出为 fields-mini.json 和 fields-macbook.json。Petersen 例精确检查了 1,022 个非平凡顶点集，十边张满无桥子图为零个，十一边张满无桥子图为 60 个，每条边在其中出现 44 次。十四位稠密槽的全部 16,384 个整数 Walsh 系数模均为 128。
