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

## 第4轮来源与核验
- 继承本仓 proofs/001-tail-and-local-rigidity.md，定位 §1–§2（变换和有限尾排除起点），基线提交3f39ed3891e587ba28c8f6add76072c1528e3f49；SHA-256: 3f099de9e6a514f855e4691a47ad2f516b79683bc77093bf7f3e897e30ceb12c. 第4轮单独给出紧性、热测度连续性和显式尾部见证证明。
- 继承本仓 proofs/003-continuous-time-realization.md，定位 §4（中间正性与原题的等价关系），基线提交3f39ed3891e587ba28c8f6add76072c1528e3f49；SHA-256: 598a303260fa4bdaa9305739889b57b11c434ea88d46f39730015d2d917e4a3d. 第4轮单独给出紧性、热测度连续性和显式尾部见证证明。
- 本轮定向检索没有提供新的可用连续刚性定理；未依赖搜索摘要作数学推论。近似数据全由本轮脚本实际生成，不把数值搜索阴性或求解器状态当作理论结果。

## 第5轮来源与核验
- 继承本仓 proofs/002-heat-positivity-and-finite-values.md（中间热质量定义与非图反例），基线41b04f504cb72fc8a87147f87e97689b205c6e92；SHA-256: 0f33ed2989a411e2bba174b96c6fdfe09273745c1b40d8edc5534a9e9e31f2e1。本轮三层构造与振荡极限另作自包含证明，未把旧核验状态当作新证据。
- 继承本仓 proofs/004-compactness-and-tail-witness.md（连续图紧性门槛），基线41b04f504cb72fc8a87147f87e97689b205c6e92；SHA-256: c360497aac046040a7e24f2dde49a7e6edceabf2fb56dec76357607e03e8d7a0。本轮三层构造与振荡极限另作自包含证明，未把旧核验状态当作新证据。
- 本轮无新增外部定理依赖；高斯积分、弱收敛与全变差收敛论证均展开。无新增文献首创性或文献穷尽主张。

## 第6轮来源与核验
- 继承 proofs/001-tail-and-local-rigidity.md §1（归一化与非负实变换），基线40f984c60677e6f4f9b8b57a06ab81002109506e；SHA-256: 3f099de9e6a514f855e4691a47ad2f516b79683bc77093bf7f3e897e30ceb12c。本轮局部积分与鞍点渐近均自包含重推，非从文献摘要导出。
- 第6轮没有新的外部定理或首创性主张。局部模型非全局解；不把数值积分符合渐近等同于原题验证。

## 第7轮来源与核验
- 继承 proofs/001-tail-and-local-rigidity.md §1–§3（归一化、尾部必要条件），基线d4504a42f1e1bf2542ea7960190114806e3d7d61；SHA-256: 3f099de9e6a514f855e4691a47ad2f516b79683bc77093bf7f3e897e30ceb12c。本轮独立推导正负贡献平衡及两项渐近，不沿用未经核验的新结论。
- 2026-11-17重新定向检索连续/光滑Cantelli问题；核对 https://arxiv.org/abs/1202.2250 的题录、v2日期和摘要。该摘要陈述的是可测反例，未提供连续版解决。未依赖搜索片段建立数学结论，也不主张检索穷尽或本轮首创性。
