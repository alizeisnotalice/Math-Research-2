# 235个指定根分配的精确排除

冻结degree54任务18..23中原先待定的235个mask全部严格产生负系数。每task54个有理根盘，共324盘；只涵盖这235个指定根分配，不涵盖其他保存样本或全部2,684,354,560个数值掩码。

生成器用100位Newton提出80位有理中心，半径10^-60；实际证明全部用整数/有理数。原证书JSON SHA `1786dc73` 前缀仅作定位，完整hash以两个独立验证输出的bindings为准；压缩文件SHA `af335a94fc729d267e904405cb236c1c53f88c2c91ac3708b32553ef349979c5`。完整Fraction验证重新计算Rouché、分组、归一化及负区间；父另独立用不同尾界及10^-100整数网格外舍入，235项全部通过，最终源码25.644秒。六个伪造覆盖/根盘/掩码/常数负项/优化调用回归全部拒绝。

父粗界证明：二元C次数D、0<r≤1时，二阶余项≤(r²/2)Σ_i i(i−1)(|Re z|+|Im z|+1)^(i−2)≤D³(|Re z|+|Im z|+1)^D r²。逐根以此上界作严格Rouché比较。系数区间每次乘法在整数网格下取floor/ceil，故保持包含；无需近似求除，只检查原始系数与非零常数项区间异号，即得常数项归一化后的系数严格为负。

解压本包证书，并解压相邻045-native54-partial/inputs-degree54.json.gz及047-native54-final6/candidate-diagnostics.json.gz到data临时目录。执行：

```
python verify_root_disks_fixed.py --certificate CERT.json --input INPUT.json --diagnostics DIAG.json --out parent-check.json
python verify_exact_root_certificate.py --certificate CERT.json --input INPUT.json --diagnostics DIAG.json --out fraction-check.json
```

不要加 `-O`；父验证器明确拒绝优化模式。源输入hash和全部任务、235mask、seed附近唯一根的分组绑定均被重新检查。全部根分配的常数项归一化积仍为C，因为原首一因子之积是C且C(0)=1。
