# 共同基域 Kummer 条件定理的有限模块检查

一般证明见 [031-kummer-two-fibers.md](../../../proofs/031-kummer-two-fibers.md)。`python3 verify_permutation_modules.py` 使用标准库精确模素数行消元，复核(n,r)=(5,3)全部81个增广向量与(5,5)全部625个增广向量；它只说明非恒定轨道张成增广模、常量线是整除情形的例外。不证明具体母根的非幂或幂类秩前提。

主线程实际运行退出0；回执SHA为10fbf0d8376fc1a23b84dbde9211af72bb5e6211d777598a4151581823a62d4b。

文献核验：2026-11-16阅读 [Kedlaya §1.2](https://kskedlaya.org/cft/sec_Kummer_theory.html) 的Remark 1.2.7及Theorem 1.2.9，确认含单位根的Kummer对应和配对非退化；另核对 [Stacks §9.24](https://stacks.math.columbia.edu/tag/09I6) 的基本Kummer前提与ζ_p扩张次数。这里的证明独立写出有限配对和置换模论证，不复制来源全文，不依赖Stacks该节其他命题。必须在共同基域K(ζ_p,ζ_q)核查非幂条件，不能直接搬用分别扩张中的条件。
