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

## 第5轮补充（2026-11-17）

- [Mande–Sanyal 2020全文](https://arxiv.org/html/2008.00266)：本轮主Agent实际阅读§3.4定理3.7、3.8及§6对应证明和相邻构造。每个谱点都有至少三对折叠的方向，以及唯一此类方向必折叠全部支撑，均属于已有定理，不能作为本轮新发现。轻锚点自包含证明用于复核；本轮正规形和通用包装推导明确以原定理3.8为入口。未声称逐项重证整篇论文。
- 本轮没有据有限计算推断一般平台谱的最优维数；此前Tarannikov的构造上端与一般上界仍严格区分。没有提交外部论文全文。

## 第6轮补充（2026-11-17）

- [Meidl 2024官方记录](https://www.aimsciences.org/article/doi/10.3934/amc.2023004)：主Agent核对题录、摘要及访问范围。摘要明确admissible分块常值空间维数等于块数，块数至多2^{n/2}，也提及非平行覆盖。Full Text入口返回登录/付费页，未取得全文。文献子Agent曾把索引片段中的Proposition 2.7及Theorem 3.4说成全文核读，来源审计后撤回；本轮不把这些条目登记为已核对全文。
- [Çeşmelioğlu–Meidl 2025](https://doi.org/10.1007/s12095-025-00805-7)，Construction and equivalence for generalized boolean functions：取得官方公开PDF；主Agent实际读取提取文本第1664–1665页Theorem 2、Remark 1及第1670页§4入口。文中partition-constant affine space仅是H，不是任意含H的共同扩展V。H的维数界不能替代本轮V的2t−1界。主证明自包含，并沿用第2轮仿射分割秩定理；未对其首次性作全球检索承诺。
- [Meidl–Pott 2019](https://doi.org/10.1007/s12095-019-00370-w)只取得题录/摘要；未取得正文，不将其原命题作为本轮已读依赖。六变量逐函数MM分类的原文访问受限，本轮不使用该分类代替共同U检验。
- [de Seguins Pazzis 2023](https://arxiv.org/abs/2307.10347)和[2026](https://arxiv.org/abs/2607.14648)的锐交替矩阵空间结论带域大小限制，不能直接用于F₂。二次子类补充改用直接ANF维数及自包含的特定pencil包含极大性证明。

实际来源文件哈希见source-files-r006.json；外部全文不归档。有限搜索阴性没有被当作一般上界；错误的四原子独立翻号推断已由24个最小码字的证书及正确Hadamard分类排除。

## 第7轮补充（2026-11-17）

- [Kolomeec 2016](https://doi.org/10.4213/mvk187)，A graph of minimal distances between bent functions，Mathematical Aspects of Cryptography 7(2), 103–110：文献子Agent实际通过MathNet原文渲染阅读8页；第105页为单个最小翻号的仿射平坦表征，第106页Lemma 1为特定MM排列交换，第107页Statement 3为相同仿射线性部的两平行陪集拼接。后两项并非一般相交翻号的E₄闭包或秩定理。主Agent核对阅读备忘录；本轮其独立重开MathNet记录页失败，未将这次失败记作主Agent阅读全文。证明自包含。
- [Kolomeec 2017](https://doi.org/10.1007/s10623-016-0306-4)，The graph of minimal distances of bent functions and its properties：只取得ePrint元数据/摘要与Springer预览，ePrint PDF当前403；未核读全文，不以摘要推断E₄或秩结论。
- [Bykov–Kolomeec 2025原文](https://math-sobolev.ru/sites/math-sobolev.ru/files/2026-02/32-3-5.pdf)，DOI 10.33048/daio.2025.32.811：文献子Agent读取公开俄文原文，其最近邻分类/计数仍围绕单个仿射平坦。主Agent未逐项重读该全文；本文定理不依赖其计数。已核查来源没有定位本文完整结构断言，不据此宣称首次性。外部PDF不提交。

## 第8轮补充（2026-11-17）

- [Filmus等，arXiv:2212.14685v4](https://arxiv.org/abs/2212.14685) Theorem 4.7 和 [Bamberg–Filmus–Ihringer–Kurz，arXiv:2211.16561v2](https://arxiv.org/abs/2211.16561) 的分割定义由文献子Agent重核：完整覆盖是原定理的假设。本轮部分族结论给自包含证明，单缺块补全后才应用完整分割界。
- [Delsarte–Goethals–MacWilliams 1970](https://doi.org/10.1016/S0019-9958(70)90214-7) 本轮仅核题录/摘要；[Leducq，arXiv:1001.2554](https://arxiv.org/html/1001.2554) Theorem 1由文献子Agent核读，主Agent也实际重开HTML核对引言、定理1及其二元特例。最小重量支撑的仿射刻画属于经典结果。归档证明另给导数与周期的二元自包含论证，不主张单缺块推论首次出现。
- [Borissov–Manev–Nikova，WCC2001原文](https://www.esat.kuleuven.be/cosic/publications/article-413.pdf) Lemma 2由文献子Agent核读：非最小码字的分解结论带有非最小性假设，不能推出全部两倍最小重量码字可拆成两个仿射块。Kasami–Tokura相关旧原文未取得全文；严格小于两倍最小重量的转述不能支付等号边界，本轮不依赖未读原命题。
- 本轮显式RM(3,6)例子仅否定“次数＋两倍最小重量”简化。主Agent完整检查其48点补集无六块精确覆盖；不把它作为实际部分分割无法补全的反例。无外部论文全文纳入提交。

## 第9轮补充（2026-11-17）

- 一般部分仿射族的法向秩≤2^c−1采用本轮自包含证明。文献子Agent检索完整仿射分割、线性partial spreads与near partitions，没有找到同一结论和归纳；检索阴性不承担首次性断言，也不把线性partial spread文献直接移用于仿射不交族。
- [Ballet–Rolland，arXiv:1203.4592v4](https://arxiv.org/html/1203.4592v4)：主Agent实际读取§4.1（HTML行319–323）和Remark 4.3的二元分支（行374–384）。第二重量3·2^(n−r−1)要求2≤r<n−1；本轮r=k,n=2k,k≥2符合。r1与r=n−1端点另有公式，未混用。[Kasami–Tokura–Azumi 1976](https://doi.org/10.1016/S0019-9958(76)90355-7)仅取得题录，未声称读其原全文。
- bent相位次数界按经典Rothaus界归属；本轮正文给基于Walsh/Poisson、任意(k+1)维坐标面奇偶性与ANF差分的直接证明，明确k≥2。参考访问到的[Rothaus原论文题录](https://doi.org/10.1016/0097-3165(76)90024-8)不计为原文证明已核读。
- 未提交外部文献全文。一般部分秩定理、近覆盖推论、有限计算与尚未解决的一般A28目标分开标记。
