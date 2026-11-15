# 第36轮主代理实际复放

已审读恢复脚本及指南，确认恢复目录须位于已挂载data盘的CodexTemp/T2下且不得预先存在。固定候选复算命令与两份算核实际argparse参数一致。

本轮在新的`/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round036-root-restore-20261115`实际执行归档脚本，退出码0，生成26项。随后独立逐项检查：解压项等于归档gzip原始字节，复制项等于归档源字节，R25重建完成流逐条等于保存results，全部SHA256一致。该重建流不声称恢复原完成顺序或原字节。详见[机器结果](../computation/results/round036/root-restore-check.json)。未重算halo或运行历史搜索。

复放命令：从课题根执行`python3 computation/tools/restore_rounds_025_029.py --archive . --output <data盘下新的临时目录>`，再依RESTORE-MANIFEST.json逐项核对归档来源。生产者报告和其独立目录历史运行保持原字节；此处单独记载主代理本次运行。
