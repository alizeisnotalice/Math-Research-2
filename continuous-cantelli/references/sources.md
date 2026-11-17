# 来源与核验范围

- 用户附件：Continuous_Cantelli_Pro_攻坚提示词_20261117.md；SHA-256: 737ebbff77478de8adbd9d996b696087f521cc2896a3ae3134ea7cf512e56b89. 定位：原题、共同起点、线 A/B。来源声称三项必要条件；本轮第1轮已独立完整推导，见 proofs/001-tail-and-local-rigidity.md。只归档自写研究材料，不公开附件全文。
- Kleptsyn–Kurtzmann, *A counterexample to the Cantelli conjecture through the Skorokhod embedding problem*, Ann. Probab. 43 (2015), 2250–2281, DOI 10.1214/14-AOP932；https://arxiv.org/pdf/1202.2250 ，版本 v2。已读电子版 pp.1–7，尤其 Theorem 1.1 后连续版说明、Proposition 2.1、§2.2 remark、式 (2.7)。论文声称可测反例及连续有界 Brownian transport 的 Root 唯一性；本次核对陈述与适用假设，未重证完整构造或 Root 定理。
- 2026-11-17 定向检索连续版后续解答，未发现可核验的解决论文；检索阴性不证明文献穷尽。
- GitHub 文件限制核对：https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github ，2026-11-17；本项目采用更小的 10 MiB 单文件归档门槛。

## 第2轮来源与证据边界
- 继承本仓第1轮，提交5c7abc333d5bd3577e1ed9ff54fe2199adb683d6，proofs/001-tail-and-local-rigidity.md 的§1、§6。原始高斯变换与有界性本轮按公式再次核对。
- 第2轮的热测度公式、有限层实现、一般正测度反例与最小噪声层定理均给出自包含证明，不依赖未核验的逆向嵌入定理。定向检索未用于宣称首创性或文献穷尽。
- 原论文 §2.2 的连续有界 transport 唯一性仍仅作为已核对陈述引用；本轮没有重证 Root/Loynes 定理，也未断言其适用前提已经满足。

## 第3轮来源与核验
- 继承本仓第2轮，提交1f8a45ca02460114764b0e01d17ddbdafbefcedb，proofs/002-heat-positivity-and-finite-values.md 的§3–§5；SHA-256: 0f33ed2989a411e2bba174b96c6fdfe09273745c1b40d8edc5534a9e9e31f2e1. 第3轮扩展有限时间实现，并实际核对时间积分公式及两层源。
- M. Beiglböck, A. M. G. Cox, M. Huesmann, *Optimal Transport and Skorokhod Embedding*, Invent. Math. 208 (2017), 327–400, DOI 10.1007/s00222-016-0692-2；https://arxiv.org/pdf/1307.3656 。2026-11-17 已读 PDF §3.2–§3.3 中 Theorem 3.8、Corollary 3.10、Lemma 3.11 及 Theorem 3.14 的适用范围。当前PDF44页，相关段落约PDF第13–16页。用于核对随机化停止时的弱闭性；本轮所需有界时间版本另给展开证明，不声称复核全文。
- 作者站点旧版 Transport2Skorokhod_4_2_15.pdf 返回404，改读上述arXiv全文；不依赖旧版检索摘要中的定理编号。
