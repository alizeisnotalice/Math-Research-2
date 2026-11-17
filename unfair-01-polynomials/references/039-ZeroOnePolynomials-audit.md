# David Zhang 的计算性声明：独立源码核验

本次取得的是公开源码，不是完整输入 data 或 proof corpus。因此，本次确认的是公开树规则和证书生成路径；不确认总次数 66 的完整覆盖，也不确认 10,501,314 个叶子或 cofactor degree <=4 的全部证书。

## 来源、版本、许可

* 原回答：https://mathoverflow.net/questions/339137/why-do-polynomials-with-coefficients-0-1-like-to-have-only-factors-with-0-1 。本次于 2026-11-17 读取。David Zhang 的回答段落中报告总次数 <=66、10,501,314 leaves、cofactor degree <=4；Zenodo <=70 是计划，不是已上传材料。此前记录的 111--141 是本次网页抓取文本的定位，不是该动态网页的持久行号。
* 源码：https://github.com/dzhang314/ZeroOnePolynomials 。
* 本次 clone：/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/source_audit_20261117/ZeroOnePolynomials 。
* 当前提交：582cdcc4a17a84531a5bda6ee5ac8d3fc0ff7766（2026-07-26T20:34:17-07:00）；完整历史已 fetch，当前没有 tag。
* 许可证：MIT。LICENSE SHA-256 为 31be45737131f82fb6a77ccb9ae4f8dcade464fc367c6bb49c80dce77da0931a。GitHub Releases 页面显示没有 release。

关键源码 SHA-256：

* solver/ZeroOneSolver.cpp：5552aede5c7aafd12f0787610d6e83bb76fcfdb6ab8b7004c7ef246362a0b417
* solver/ZeroOneSolver.hpp：c21cd56b25da91c0f97112ecde716bd8e8ef0486aa35968ed4f05291bea8185f
* solver/PositivstellensatzCertificates.jl：73e6c4ed87b13f1fb90443e8aa2f6696d0bf05dca8c57d5b5db13292d60d8bd0
* solver/DistributedVerifier.jl：b817cb0b201dcdae7ef7c205e1ce35303dabb2cc2f19a279e58e0a7c45d0ad27
* solver/data.sha512：67a06a65ad3add08970fb513fa2727d2cbbdc1c6d9f8215618e8cc9f5ba7be1f

## 原声称与本轮验证

| 项目 | 原作者声称 | 本轮实际验证 | 结论 |
| --- | --- | --- | --- |
| 总次数 <=66 | MO 回答报告已覆盖 | 只有源码和外置 data 哈希清单 | 未独立复验 |
| 10,501,314 leaves | MO 回答报告 | 没有 data/leaf/proof corpus | 未取得 |
| degree <=4 证书 | MO 回答报告 | 读到生成代码，未取得证书 | 仅核验实现路径 |
| 浮点 LP 后精确核验 | MO 回答说明 | 读到 HiGHS -> FLINT fmpq -> 精确稀疏矩阵重乘代码 | 源码层面确认 |
| Zenodo <=70 | MO 中未来计划 | Releases 页面无 release，仓库无 Zenodo 链接 | 未取得 |

## 树规则和精确证书接口

ZeroOneSolver.cpp 第 495--524 行从 x^M、x^N 系数产生每个 p_i 的零/严格正 case split。第 366--484 行依次处理已知二元变量、零 RHS 双线性项、单项 0-or-1、domination 和最短未定 RHS 方程。ZeroOneSolver.hpp 第 289--367、458--470、527--558 行维护严格正性并使用非负和为零的简化。ZeroOneSolver.cpp 第 513--535 行只在无规则可用且仍有未知变量时输出 leaf；反向 case ID 去同构由 exceeds_reverse 完成。

DistributedVerifier.jl 第 32--45 行按 d=1,2,3,4 构造 LP。PositivstellensatzCertificates.jl 第 368--408 行调用 FLINT 的 fmpq_mat_can_solve_fmpz_mat_multi_mod，拒绝负系数，再用 Julia 的精确有理稀疏矩阵重乘并要求单位向量。该代码确实实现浮点支持定位、精确有理解、精确再相乘。

限制也很明确：CheckProofDirectory.py 第 31--58 行只检查目录 shard/文件名；DistributedVerifier.jl 第 59--69 行对已有 proof 仅比较其开头的输入字节并跳过重算。当前仓库没有独立 proof-text parser 来重验已写出的有理等式。

## 可得材料与覆盖清单

.gitignore 忽略 solver/data/ZeroOneEquations 和 solver/proofs；本次 clone 中两个目录均不存在。实际运行 CheckProofDirectory.py 退出 1，报 FileNotFoundError: proofs。Julia 未安装，因此没有运行 DistributedVerifier 或任何证书。

data.sha512 有 1120 条外置数据哈希，不含 data 本体。按 DistributedVerifier.jl 第 49--50 行的任务集（m=5 或 m>=7，m<n）复算：

* 截至 45：291 个期望对、291 条目；
* 截至 66：758 个期望对、758 条目；
* 截至 70：872 个期望对、872 条目；
* 截至 78：期望 1124 对、实际 1120 条，未完整。

这只证明哈希清单的文件名覆盖了该程序自定义的任务集。ParallelSolver.py 第 49--50 行称 m=1,2,3,4,6 已解决，等次数也未进入该枚举；相应证明材料本轮没有取得。README 当前仍写 <=45，故不得把源码或哈希清单升级为本项目的 <=66 或 <=70 有限排除。

## 最小实际运行

使用 Apple clang 21、未改源码，已实际编译并运行两个全 case-ID 小例。以下 stdout 和二进制哈希明确为不同文件：

* M=5,N=7：退出 0。stdout 路径为
  repro_build/tree_5_7.stdout，0 bytes，SHA-256
  e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855。
  二进制路径为 repro_build/tree_5_7，56,248 bytes，SHA-256
  fd52b42d7c69cebe84fbe3ee5aa3b48bfe75bb87ac51d7d0b6a5e8d5f1d4cf2c。
* M=7,N=13：退出 0。stdout 路径为
  repro_build/tree_7_13.stdout，0 bytes，SHA-256
  e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855。
  二进制路径为 repro_build/tree_7_13，56,456 bytes，SHA-256
  4579f6dda6129183c80396c1f8859895a6cafef9ddfa973774cbed67e865131a。

这只说明当前树程序在两个小参数下可编译运行，绝不构成对 <=45、<=66 或 <=70 的独立复现。最小复现脚本见同目录 repro_minimal.sh。若将来取得 data/proofs，应先执行 sha512sum -c solver/data.sha512，再实现或运行 proof 文本的精确等式重验。
