# 文献核验记录

日期：2026-11-16，Asia/Shanghai。实际工具：web 搜索、打开原文、页内定位。未声称检索穷尽。

查询：`Fourier dimension sparsity sqrt conjecture Sanyal improved bound 2025 2026`；`"Fourier dimension" "sparsity" "2025"`；对应2026查询；`Khalyavin Lobanov Tarannikov Fourier dimension sparsity`。

| 来源 | 本次阅读范围 | 核验结果 |
| --- | --- | --- |
| [Filmus等题集](https://simons.berkeley.edu/sites/default/files/openprobsmerged.pdf) | 定义及PDF第10页L0问题 | 一阶系数全部非零；{0,1}编码；John Wright；不是仅要求依赖全部变量 |
| [Sanyal 2019](https://theoryofcomputing.org/articles/v015a011/v015a011.pdf) | 定义、定理1.2、引理1.3、证明思路；第8–10页陪集与求和论证 | 普遍界含log s；核对共同查询空间、分支谱并集与调和和；未逐项重证被引用的单色仿射子空间引理 |
| [Mande–Sanyal 2020](https://arxiv.org/abs/2008.00266) | 题录与摘要 | 自适应奇偶树、折叠性质；不能把树深界改称谱维数界 |
| [Khalyavin等2016](https://www.mathnet.ru/eng/semr755) | 题录和摘要；PDF链接部分请求跳回记录页 | 平台函数与同谱支撑构造。未凭摘要认定本轮树构造新颖，未引用未核验的一般秩定理 |
| [Ghosh等2026](https://proceedings.iclr.cc/paper_files/paper/2026/file/634cb3ace86e0908a721865a73f2e36e-Paper-Conference.pdf) | 题录、摘要及PDF第4、6页引用Sanyal处 | 两处写O(sqrt s)，所引Sanyal定理实际为O(sqrt s log s)。此引用不提供去对数证明；本次不评价文章其他定理 |

上述条目身份核验与数学证明核验分开。检索未提供可核验的一般无对数证明或无穷反例族，不据此宣称世界范围内没有新结果。

归档大小规则实查：[GitHub官方说明](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)，普通Git超过50 MiB警告，超过100 MiB阻止；本课题自行保持逐文件小于10 MiB。论文全文不提交，只有自写论证、题录和必要精确证书。输入附件的来源哈希见 input-source.md。

## 第2轮补充

- [Irreducible subcube partitions v4](https://arxiv.org/html/2212.14685v4)：2022初稿、2026-09-08版本；作者表与版本在arXiv摘要页核验。实际阅读§4.2.1定理4.7完整证明、命题3.12及紧性定义。结论适用于紧仿射分割，不要求不可约。正文给出自包含线性代表证明与本题应用。
- [Affine vector space partitions v2](https://arxiv.org/html/2211.16561v2)：核验题录、紧性框架和二元下界入口；主证明采用上条可定位定理。
- [Tarannikov 2006记录](https://www.mathnet.ru/eng/dm65)：主Agent核验题录摘要，构造子Agent核验原文定理2、3及后续猜想。构造某个秩区间不等于证明一般最大秩；本报告不使用其未核验的一般锐上界，不排除全部平台函数类。

子Agent报告均经主Agent数学整合；原题的完成状态没有因子Agent的“审计通过”而改变。没有用后续论文的误引替代原定理。

## 第3轮补充

- [Tight Chang’s-Lemma-Type Bounds for Boolean Functions，arXiv v2](https://arxiv.org/html/2012.02335)：主Agent实际核对权重依赖界与§6.3组合引理6.8；[FSTTCS 2021正式记录](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.FSTTCS.2021.10)用于作者、题录、DOI及会议版定位。一般权重界仍带log s。组合稀疏度式1+t²(s−1)的原假设蕴含非零均值，不能用于平衡字符；本轮证明显式用δ₀=1[ĝ(0)≠0]修正并核验两类种子。正式PDF请求超时，未声称通过该请求读取PDF全文。
- 文献子Agent检索共同查询、权重预算及寻址组合，未定位到与本轮任意重复叶字符树定理完全相同的已核验陈述；此阴性检索不构成新颖性证明。文献子Agent初稿曾遗漏组合引理的非零均值条件，主Agent以平衡字符发现并修正。

## 第4轮补充（研究2026-11-16，归档2026-11-17）

- [Khalyavin–Lobanov–Tarannikov 2016原文](https://m.mathnet.ru/php/getFT.phtml?jrnid=semr&option_lang=eng&paperid=755&what=fullt)：本轮成功取得23页PDF。主Agent核对印刷页1348的Titsworth卷积判据、1355的递归族平移轨道等式、1357–1358的命题4/定理6及1360–1363的一般块谱讨论。此前轮次只读摘要的范围不沿用为本轮全文已读；本轮也未逐项重证文中全部计数。地址合法符号的“置换×独立符号”分类属于已有结果，本轮证明其共同线性符号空间的秩限制和Sidon等号频带闭包。
- [Tarannikov 2006](https://www.mathnet.ru/eng/dm65)：文献子Agent核读俄文全文，主Agent核对题录及提取文本定理2、3与最后猜想。定理2是仿射秩的构造区间2h至2^{h+1}−2；一般上界是定理3的2^{2h−1}−2^{h−1}+h，不能把构造上端当一般上界。支撑16的仿射秩4、5、6属于已证明的特殊情形。
- 搜索方向包含同支撑平台函数、谱符号切换、仿射bent族及向量值bent构造。Hodžić等[arXiv:1811.04171](https://arxiv.org/html/1811.04171)由文献子Agent核读定理3.3、4.3与相邻开放问题：可作合法符号/非仿射谱构造工具，未获得超过平方根尺度的谱秩证据。本轮主证明不依赖该文。
- 下载来源哈希见`source-files-r004.json`，只登记可读PDF，不公开外部全文。最初名为`semr755.pdf`的本地下载实际是HTML响应，已检查文件头并排除；最终引用成功的`semr755-download.pdf`。下载成功和哈希一致仅证明来源文件一致，不证明数学结论。
