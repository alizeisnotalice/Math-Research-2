# 来源与继承状态

准确来源文件、实际本地提交、SHA256、阅读范围见 [sources.json](sources.json)。旧专用克隆实际位于提交 `9079302d4f9fb8b75630f2405466352dd72104c8`；它与附件目录所列总仓库基准不是同一提交，不能只按时间合并。

- `legacy_input`：51原子准确有理配置，逐字节复制。JSON 中历史 `float_ratio` 是旧显示元数据，任何本轮数学计算均不读取它。
- `legacy_interval`：旧第4轮上下界分数，本轮核对算术；本轮没有重跑旧核上界。上界保持“原来源已证、本次未重算”的状态。
- `atomic_foundation`：有限原子矩形公式及固定N达到性相关章节已阅读；第1轮周期认证使用本轮自足证明。
- `legacy_area_algorithm`：主算核采用旧整数前缀质量与扫描线方法，新增全半径尾界、完整窗口、外包围和逐项根号认证。独立算核另用逐列累加质量和条带并区间方法。
- `reflection_draft`：旧第5轮反射草稿已独立文本审读；它未包含在当前专用克隆第4轮提交中，单独绑定文件哈希。历史实验未写成本轮实验。
- `uniform_deletion_counterexample`：反例排除“任意指定原子均可用统一 cε 控制删除损失”，不排除带额外选择规则或尺度状态的命题。

公开文献入口按作者原文核对：[Tomar](https://www.mefyi.com/papers/weak-type-squares/)、[Lin](https://coolrmal.github.io/articles/two-partial-balayage-principles/)、[Melas 2003](https://annals.math.princeton.edu/2003/157-2/p08)。引用精确命题时须另绑定实际版本和阅读页码；网页入口存在不等于全文已核。未将第三方论文全文纳入本仓库。
