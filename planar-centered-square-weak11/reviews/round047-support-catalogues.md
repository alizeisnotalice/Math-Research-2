# 第47轮候选完整性、外竖支持与上游域独审

结论：PASS（本报告分工范围）。数学稿SHA256 `3e391beace80ec39a668eaeb81de218d8ab2d982da39d919f8f3ead749d5a3b7`，新增模块 `overlap_position_radius.py` SHA `929cf92baca532b92159ff81542390acbff3688386dd5104decd0ee55a8ed524`。独立审读固定完整candidate、外y支持/同y事件/引用计数/Fx常驻和上游有效域前提；相邻异速gap证明、min树及β/旧δ的独立重放由envelope另审合并。本报告不重复旧冻结枚举器全文/旧小例，不运行44M/45M候选bbox主循环。

## 数学及代码的完整性边界

固定y区间端点把F_y分成正高条带。在一条带共同出现的盒，真实竖支持有正长交集；对这些盒和F_x边界保持横端点全预序，足以固定该条带的一维并区间结构，从而面积仿射。候选自身的两端点与F均参与，因此当前横向为空的候选也被约束保持空，不能先取新中心活跃子集。无正高y支持的候选可以永久删除，端点零测集不增加面积。

源码对输入全部boxes逐项建family，不按横长度或当前clip活性过滤；重中心CLI明确取candidate_boxes，未取新活盒子集。只删除外y支持已证非正的盒，且记录索引/拒绝原因；完整x节点仍含这些raw端点，中心排序要求因而可能更保守，不造成漏约束。固定y速度被检查为0，全部横速度须为0或1。

对y0用定向下界、y1用定向上界后再与F_y裁剪，所得外区间包含真实y支持；真实正长相交必在某个正高外条带共同激活。仅接触的外区间没有正高共同条带，相同y处的全部进入/退出先完成，再读取下一正高条带。每个支持盒的两个canonical横节点都引用计数，允许多盒共享以及同盒两端重合；一个盒退出不删除仍被其他盒引用的节点。F_x左右节点始终另有1个基础引用，即使与盒端点共享也不被退出事件撤掉。

上游domain、H及原来未列入候选目录的空盒余量是不可放松前提。原slab输入直接使用H和delta_empty；recenter输入只使用remaining_source_H和remaining_original_empty_radius，并保留全部原候选。新模块不自动证明这些remaining字段的来源，调用者仍须独立建立原source/omitted安全域减去重中心距离后的余量。旧global raw半径只用于诊断/比较，未重新成为final硬帽；final明确取remaining source、remaining original empty与新overlap半径的最小值。输出保留external_premise，不能单凭status把未经验证的原目录或重中心域认证为完整输入。

## 本次实际验证与五个新slab

本次独立脚本实际9.390542秒，CPython3.15.0，单进程、exit_code=0，无随机种子/浮点容差。六例索引SHA `185064c894cc24fb18db72da588b1ab040d50885dc54e3ec1f096663f780d4a3`；五个新中心索引SHA `f6863ff64889f91b261cff366af33b72927e95505541620c7bfc9be123a406b2`。R39原slab按已验hash `c0d262c39a01db8d50cd2fee1a00c7867aa8e08417c4fc85522103b70ce13b25`复用已审目录，不重复旧实例计算。

5个bracket-refine虽然暂存在round045路径，均为本轮新增实例，明确不在已封第45轮11例中。本轮复用冻结独立检查函数 `review_round045_enumerators.py`（SHA `5075bf81a2c0c0c137e1c305395c2254f047943c415ee0d495f84bd6353a2279`）实际核验：R41基点SHA `ab714c7f457ccf69f68fd629c947021faab16631c8995caa16b67c1ca2752db3`、原26/27标签对应五方向速度、物理t=P8·s和H=P8/100000000，其中P8=`328734960629/100000000000`。新中心归一化参数q0/q2/q3为−139/172800000，q1为−21/57600000，all为−45/57600000。5份manifest及对应历史runner版本均匹配。

全部5份3953副本目录、source/span/reach事件、候选对预筛选、访问数恒等式、133259活盒质量/半径/四个原始端点及拒绝摘要见证均准确通过；实际H=delta_empty=P8/100000000，统一R=10普通尾界<1。候选访问计数合计232049500，包含活盒；未逐盒重放该大分类循环，省略盒的穷尽分类仍依赖已审冻结枚举器、其精确判据和本次完整源/候选对检查。

## 六个外支持目录及引用计数

独立用整数平方和dyadic定向根界逐项重建六例全部外y支持与恒空记录，不调用生产min树；原candidate数组哈希、全部canonical横节点集合、相同y事件分组、正高条带数及所有引用收支均吻合。每个正高条带F_x的两个基础节点始终存在，末态只有F_x基础引用；真实数据canonical节点共享引用最多65，故引用计数确有实质作用。报告限制pair在其共享外条带有有效活跃owner，给出的gap由两个端点的定向界准确得到。

| 例 | 完整candidate/正外支持 | 正高y条带 | final扩大倍数 | final限制来源 |
| --- | ---: | ---: | --- | --- |
| R39 | 25466 | 9319 | 约5.98 | overlap端点 |
| q0 | 26331 | 9670 | 约4.98 | overlap端点 |
| q1 | 26331 | 9670 | 1 | source H |
| q2 | 27133 | 9998 | 1 | overlap端点 |
| q3 | 27133 | 9693 | 1 | overlap端点 |
| all | 26331 | 9670 | 1 | source H |

合计158725个candidate全部在外支持目录中，无实际Y恒空删除。旧/新半径、final三项min及准确倍数与六例索引逐分数一致；只扩大R39和q0两例的final，不能写成六例全部扩大。这里的半径仍是充分范围，不等于最近有贡献事件或面积饱和墙；外支持还可能增加保守相关边。

另实际调用新模块作5个有界行为例，并以独立外支持/事件引用检查复核：当前x为空候选仍保留并受自身未来激活比较限制；共享端点在一个盒退出后继续存在；与F_x共享的节点保持基础引用；传入remaining H=1/40、original empty=1/200时final不超1/200；固定Y恒空正确记账。全部通过。这些人工例不等于已验证某个真实recenter文件的剩余域来源；本轮六真实输入均为原始slab schema。

## 文件、复用与尚未完成范围

独审脚本 `/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-r001-audit/planar-centered-square-weak11/computation/review_round047_support_catalogues.py`，SHA `29d29f5de3d28309c681d5e267b2af13cb2c71512e2038a7641074f2106818e5`；实际输出 `/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round047/independent-support-catalogue-review.json`，SHA `9c9df22d3ba061a21a47ba0ba40464a4a5b8d6df18c02e4400370886e09c1735`，保存六输入/输出hash、五新目录详情、支持/引用计数、全部半径及实际命令/版本。标准库依赖，另复用上列已归档R45独立checker；五个人工小例动态导入已绑定的新模块，不对真实6例重新运行生产树。

六例最小mixed gap、树合并与真实β/旧δ由envelope本轮独立list/heap及sort-merge复核合并。本组件结论不替代那部分证明或计算；新模块不产生一般C₂上界，也不单独给新的下界候选或多维最优性。未修改作者代码或main/T，未提交/推送。
