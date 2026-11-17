# 来源、改动及证据范围

本轮公开编号35。工作目录名 r037–r041 为并行临时子任务号，不代表已经公开的轮次。
本包仅收录本课题自写证明、代码及最小有限结果；不含文献全文、凭据或远程连接配置。

- gap-mixture 来自工作子任务 a022-global-gap-mixture-r037；全n输入B/A1发散和B_q上界
  分别已在公开第34、33轮归档，纤维交错式见此前第22–23轮。
- weighted-reinsertion 来自 a022-weighted-reinsertion-analytic-r038。公开副本修正一处LaTeX
  `qquad`漏反斜杠；分析脚本原样复制，实际主线重放另记Python3.15.0。
- completion-transfer 来自 a022-weighted-average-r039 的修订版：删除未获支持的早秩
  质量实现说法，保留精确rho秩比较因子及明确的未证质量分布。没有复用该旧说法。
- comb-obstruction 来自 a022-weighted-interval-geometry-r039。公开版明确n≥5、修正LaTeX，
  移除未用于证明的非公开小调度程序引用；结论只覆盖规定的早comb骨架顺序。
- joint-geometry 来自 a022-joint-geometry-kernel-r041。公开版将无穷延拓端点明确写为
  (θ,u)=(0,0),(1,0)，内部u=0仍取有限值；其余数学内容及审计脚本不变。
- spatial 保留 a022-n7-spatial-crt-r035 最终选定目录的原嵌套结构、源代码及内部清单。
  n7原始40模文件是已完成生产的记录；本轮重放读取它们，不重复完整生产。
  原distance_counts的Manhattan语义由附带r018原始C++源确认；平方欧氏总和另存字段。
- small-gate既保存准备子任务实际产生的60模，也由主线在最终树中重新编译并运行60模，
  以独立完整联合表核验25字段。其结果不能冒充全n7生产的第二次独立枚举。

主线复现结果见final-tree-replay-receipt.json和parent-small-gate-exact.json；每个
子进程均实际wait/reap并确认PID退出。临时完整stdout/产物不批量进入公开仓库。
在途CPU/GPU生产按独立任务账本维护，不纳入本包已完成结果。
