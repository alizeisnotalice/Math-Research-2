# 六点总次数50完整认证

50层68,405分片全部独立语义重算：1,724,108,390原始支撑对、819,019,379规范代表，其中814,369,357排除、4,650,022公平、未知0、错误0。结合此前≤49，严格有限范围推进至每侧≤6点且总次数≤50；每侧≤9点仍≤34。

父check_total50.py另独立重建49/50输入并集123,859片和旧58,036／补65,823的无重无漏分割。完整读取本层44,976,106,172解压字节、3,134,349,873压缩字节，逐片核输入、内容hash、行数和语义报告绑定，实际140.3305秒通过。这是完整内容与覆盖核验，不冒充另一遍语义重算。

实际Python3.13.18。run50.py冻结源9e04a011fabac2f66dd9942178a3af91b79f515121e392e77711dddca6ef406f；实际独立语义验证器a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b。完整报告无损压缩，原SHA82641ba117b9c5dcb35abe914041692491649b0d52eda3536c62d384b2c436aa。observed-running50.json保存启动冻结绑定，完成时再次核验。所有本层验证进程已退出。

## 复现与共用材料

输入计划、旧补events、源程序、start-binding49.json、manifest-bound.json和check_union.py复用相邻[045-support49](../045-support49/README.md)，不复制同一套文件。将其中.gz解压到工作目录，按计划重新生成原始分片账本，再运行本包run50.py（其冻结本机路径需映射至对应工作目录）及check_total50.py --binding <start-binding49.json> --report <解压后的semantic-report.json> --out <输出目录>。后者observed路径也需映射到本包observed-running50.json；路径映射版本应另存，保留冻结来源和原hash。audited_chunk_verifier的verify_chunk接口用于本批，顶层旧45/46入口不用。

原始可再生账本位于CodexTemp的round041-macbook-support50/run及round041-macbook-support50-resume/run，本层约3.13GB压缩不入Git。本包保存逐片绑定、全部语义结果与复现源，不把仅读取passed字段当作验证。51及更高层、本题任意次数均不由本包结论推出。
