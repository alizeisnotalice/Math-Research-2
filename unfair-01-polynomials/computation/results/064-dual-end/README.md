# 双端分歧定理及七个二重余根

[任意次数定理](../../../proofs/064-dual-end-ramification.md)给出两端首个分数偏移 \(r,s\) 的大素数分歧条件。[恰七次的边界](../../../proofs/064-seven-double-residues.md)给出两个七次式参数族和模 \(p^2\) 条件。[独立审查](independent-audit.md)核对局部子域、对数尾项和合流矩阵，特别核对了在局部域中求解后再取归约。

精确算术使用 Python 3.13.18 和 SymPy 1.14.0。

    python3 -B verify_seven_residues.py --out /Volumes/data/CodexTemp/T2/a30-seven-residues-new.json

该运行覆盖 162 和 1,277 个平方自由且常数项非零的七次式，2,878 个整数多项式合成算术样本通过。这些样本用于验证模提升运算。

对一个平方自由的次数 64 二元积，若 \(p=163\) 的 gcd 次数恰为七次，第一道检验要求负倒根多项式模 163 被参数族中的 \(J^2\) 整除。C++ 分开枚举二进制系数的两半，以整系数余式为键，无近似运算。脚本对小范围模型与 Python 直接多项式除法逐项比对，试跑覆盖正命中和空结果案例。MacBook 上次数 64 全参数任务已启动，完整覆盖和命中的模 \(p^2\) 验证仍在进行。

    clang++ -std=c++20 -O3 search_j2_degree64.cpp -o /Volumes/data/CodexTemp/T2/a30-search-j2
    python3 -B verify_j2_smoke.py /Volumes/data/CodexTemp/T2/a30-search-j2
