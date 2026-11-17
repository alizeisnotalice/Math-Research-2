# 第32轮来源与核验范围

- Sepehr Nezami、Jeongwan Haah，*Classification of Small Triorthogonal Codes*，arXiv:2107.09684v2 (2022)。[原始PDF](https://arxiv.org/pdf/2107.09684)。Corollary 6（PDF第8–9页）、Section IV的完整分类说明、Table II（PDF第18–19页）。仅导入weight32且affine span≤8的九个代表，表中6–14行；分类搜索依赖原文，本轮未重跑。主Agent直接打开原文核对范围，两个agent独立核对代表及低秩压缩。
- 缓存PDF SHA256：`b511847558fd583ec422ef6e2a330457b00dd833996676e06c746da044c2765d`。缓存文件旧名字含bravyi-haah，但实际内容和作者是Nezami–Haah；未将PDF归档入Git。
- [作者代码](https://github.com/sgnez/Tri_from_RM/tree/89c7c5fc8b871eec1eddda10c4043106c553de45)，MIT，commit `89c7c5fc8b871eec1eddda10c4043106c553de45`，`Code_Builder.ipynb` execution_count319。最小输入 `computation/round032/rm48_weight32_representatives.json` 保留各代表的表号、原始行定位、变量编号及嵌入约定；仅抽取数学公式，不导入完整源库。
- 第19轮六变量全投影证书作为八变量剪切推论的历史依赖；本轮运行记录给继承文件哈希，不把哈希核对说成重跑证书。

新的核性质、无穷外包族及signed枚举结论由本项目证明/程序承担，并非原文声称解决A28。常用bent次数界与平衡导数性质沿用前轮固定定义；本轮两个种子的全部非零导数实际检查平衡。
