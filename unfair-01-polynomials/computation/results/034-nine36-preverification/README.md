# 44／36运行层中已返回的95个局部支撑证书

本包不声称44／36层或任何分片的完整覆盖已经核验。收集时只有nine36-m15、m16、m17三个搜索完成，各自的独立逐项验证仍在运行；其它37片尚未开始。

这三个搜索返回的95个支撑均已取得精确局部证书：92个严格对角注入、3个正多项式恒等式。主线程独立重新检查全部证书，并将每个支撑及行号与本包的原始unknown快照逐条绑定；三条正恒等式由不导入搜索器的Fraction验证器从支撑重建。该结论对这95个具体支撑成立，不依赖搜索有没有覆盖其它支撑。

发现器的历史字段search_summary_sha256命名不准确：它实际保存summary中的未压缩账本哈希，并非summary文件本身哈希。主线程最初按文件哈希比对时拒绝，查明来源代码后在最终验证器中分别绑定文件SHA与声明的账本SHA，未篡改旧输入，也未将未执行的全账本验证说成完成。

## 复现验证

    python3 verify_95_supports.py --topic /Volumes/data/CodexProjects/T2/unfair-01-polynomials --bundle . --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/nine36-95-recheck.json

line15/16由既有positive_certificate_search.py取得：乘子次数1、正原子次数2、最大有理重构分母4096。支撑分别为：

- I=0,1,3,4,7,10,16；J=0,1,2,5,7,11,14,20。
- 同一I；J=0,1,2,5,7,8,11,14,20。

实际环境为Python3.13.18、NumPy2.5.4、SciPy1.18.1、SymPy1.14.0。原输出保存在positive-preverification目录；发现器与验证器分别归档于此前轮次和computation目录。

line8的支撑为I=0,1,2,3,5,6,9,16、J=0,1,2,4,6,8,9,12,19。使用第33轮稀疏构建器、乘子/正原子次数均2，得到2381×5168整数矩阵。LP选282个活动列后精确有理求解，重算全部整数等式并检查非负性；最后独立Fraction恒等式验收。整数矩阵SHA为c0167e38d2b17599f3294e871526b4860baf5b43d8a469488727c2c96b6d5af7。矩阵可再生，未将大数组写入Git。

重现发现过程时，将本包recover_line8_cached_active.py复制到新的data临时目录。在该目录下使用第33轮matrix-builder/build_matrix_cache.py、以上I/J与degree2/atom2参数生成line8-d2a2-matrixcache子目录，再执行恢复脚本；它拒绝覆盖已有输出，自动调用独立验证器。参数、输出格式与精确活动基均保留。

主线程另对line8正确证书作生产回归：接受原证书，拒绝10^-30级精确扰动、负正原子权重和漏严格变量。对应测试程序沿用第33轮test_last_certificate.py。发现和局部验证进程已结束，父批次的三个完整账本验证器继续运行。
