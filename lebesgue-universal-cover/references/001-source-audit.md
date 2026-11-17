# 第1轮来源审计（2026-11-17）

## 已读取入口与证据范围

- Baez–Bagdasaryan–Gibbs，2015：[原文](https://arxiv.org/html/1502.01251)，§1、§3。模型允许反射，Reduction 3 实际使用反射。文献声明上界0.8441153。本轮审查定义与相关论证段落，没有重算完整构造。
- Gibbs，2018：[An Upper Bound for Lebesgue's Covering Problem](https://arxiv.org/abs/1810.10089)，摘要声明上界0.8440935944。它是提示词更强上界的直接来源；完整万能性证明及面积认证尚待逐步复核。
- Tomar，2026：[作者入口](https://www.mefyi.com/papers/lebesgue-universal-cover/)，论文p2/Thm1.1；[程序仓库v9](https://github.com/tejasvi/lebesgue-cover-bound/tree/v9)。来源声明五体D、B3、B5、B7、B9下界0.8366。模型含所有欧氏等距；这些测试体各自反射对称，因此本下界也适用于仅正向等距版本，但不能据此断言两种万能覆盖常数相等。
- Zeng，2026：[原文](https://arxiv.org/html/2609.01284v1)，S2、S4–S6、S14。层级包含连续有限维的Reuleaux形状族，非有限离散测试族；来源声明单调收敛及O(M^-2)误差。未逐条独立重证其定量逼近定理；不依赖它证明本轮3/n有限网。平面定宽补全在S2援引经典凸几何文献，尚未逐本核验；本轮证明不使用该约化。

## 下界包实际检查

克隆v9并读取程序/证书说明；commit为 `10af98606c98ca9ffd0d9f402cae1e2efe74b792`。
README.md SHA-256为 `3ddb15d188553763cc108fb7c7a4a97d2838800e88fd081845b5c9345d6fbcf5`。
README行21给五体局部放置浮点面积0.8372950064；没有定向认证，不能当作严格上界。
README行61给11D全部parts球算术检查估计约395小时；行154–155给完整Lean内核
与native构建的巨大资源量。本轮没有运行这些完整验证，也未下载大体积补充证书。
默认本机Python没有所需NumPy/flint等库；不为初始化盲目安装或重复长期计算。

下一步可先审计 `independent/opt11_point.json`，对明确放置建立向外包络与严格面积上界，
判断现有五体族的能力上限。即使证明该族最小值小于某数，也不能推出Λ小于该数。

## 发布与来源权限

仅归档自写摘要、脚本、结果及来源定位，不公开文献全文或第三方证书全集。
GitHub文件限制核验：[官方文档](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)。
项目仍采用单文件10 MiB默认上限；本轮文件均远小于该值。
