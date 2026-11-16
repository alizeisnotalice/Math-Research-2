# H10(3,5) 精确认证包

本包认证 H=1+x+x²+x⁹+x¹⁰、C=H(x³)H(x⁵) 的全部归一化非负实分解。完整必要布局为何恰17类见 [证明](../../../proofs/031-H1035-theorem.md)；程序验证有限域前提、母盘、完整布局及严格系数符号，不把数值根近似当证书。

实际 Python 3.13.18。验收只需标准库；上游素数输入搜索使用 SymPy 1.14.0。本轮独立验收的是选定的两份正证书，并未独立重算所有不超过5000素数的阴性分类。`H10-double-cycles.json` 记录上游完整搜索输入哈希。

在本目录运行；输出回执请放 data 临时目录：

```sh
python3 verify_h10_modular.py H10-mother-certificate.json H10-double-cycles.json
python3 search_h1035.py H10-mother-balls.json /Volumes/data/CodexTemp/T2/unfair-01-polynomials/h1035-search-recheck.json
python3 verify_h1035_coefficients.py --balls H10-mother-balls.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/h1035-exact-recheck.json
python3 verify_h1035_layout.py
python3 test_h10_modular.py
```

后两项在自身目录更新同名回执，可复制整个包到临时目录复现。`verify_h1035_layout.py` 检查归档的搜索和独立验证结果；若要对照新生成结果，先在临时副本中将它们分别保存为 `H1035-search-witnesses.json` 和 `parent-coefficient-verification.json`。

`extract_h10_mother_balls.py --source ../026-H10-exact/H10-a2b3-root-balls.json --output <新文件>` 只提取五个母盘，不采用旧2/3提升根。母盘验证器从固定H的系数重新计算有理Rouché界并验证共轭十盘互斥。母根隔离继承输入、重新验证；历史来源标签不是证明。

搜索器采用96位定点向外取整；独立验证器采用Fraction区间，全量重建17类。15类严格负系数，余恒定与标准类由精确整数恒等式确认。`negative-witness-summary.json` 给出每类一个向上放宽的百万分母见证，最大上界−8869/12500。区间跨零项不按零处理；两合法类的非负性由整数表达式承担。

搜索四位mask的1表示5纤维，独立验证器四位pair_3_bits的1表示3纤维，两者互补；标准分别为0与15。布局验证器从全部1024个原始比特组合重建34个平衡向量和17个交换类，未信任成功标记。母盘、系数、布局和模分解语义篡改回归均已实际通过。源码与实际回执哈希见 `parent-validation-receipt.json`。
