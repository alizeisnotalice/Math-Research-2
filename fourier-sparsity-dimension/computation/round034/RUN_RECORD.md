2026-11-18，Python 整数 Walsh 和 GF(2) 消元。归档源码在本机实际运行。复现时在本目录依次执行：

    python3 verify_amplitude.py ../round032/exceptional_signed_mask_result.json amplitude.json
    python3 audit_seed.py ../round033/mm_completion_result.json seed.json
    python3 evaluation_ranks.py
    python3 kernel_profiles.py

前两项核验全部19个局部三元组、64个图饱和符号、256个相位与显式十维共同码的1024个相位。消元程序重建九类支撑，逐一计算限制单项式的秩。12进程精确核占用枚举共检查97155个子空间，192个使指定对偶块投影单射。本机一次 top 采样整机0% idle，枚举结束后进程退出。

MacBook 的 Swift 源与 Metal 结果已归档。GPU与CPU直接字符求和比较256×128个系数，逐项相同。只用整数运算。解析无穷维结论由证明承担。早期临时Wolfram可行性模型限时结束，没有提供排除结论。

MacBook 随后以15个进程枚举七维空间的仿射平坦。四维方向11811个，陪集94488个，S32与四维平坦最大交数12。此次任务CPU累计4.096秒，墙时1.108秒，平均3.696核等效。采样窗口整机最低0% idle。15个worker均退出。GPU批次与这组CPU全核作业分别计量。源码 affine_flat_audit.py、R34S32AffineMetal.swift 及两份远端JSON保留了实际运行结果。
