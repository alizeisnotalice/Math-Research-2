# 第35轮精确核验，2026-11-18

在本目录实际运行以下程序（Python 3.15.0，整数，未用随机数或浮点容差）。

    python3 exhaust_4d_local.py four_variable_local_result.json
    python3 verify_affine_double.py affine-double-audit.json
    python3 all_bases.py

第一项完整枚举固定基相位的390个合法八点块，并保留两种取等见证。第二项固定超平面，遍历896个基相位和各256个切割。第三项以12进程扫描全部896个四变量bent基相位，每一合法八点块均检查全部70个半块；累计26880个仿射块实例各有6个合法半块，322560个非仿射块实例各有2个。两种类型另含空翻号与整块翻号。

本机 Apple M6 及 MacBook Apple M5 Pro 的 CPU 与 Metal 另行比较六变量非仿射和仿射16点块的全部65536种局部修改，各4194304个系数，得到4和8个合法翻号，并重验八变量种子的512个谱系数。整数结果均零分歧。八变量单块维数1使用第34轮已核验的无四维平坦证书及本轮解析分类。全维推论由解析证明承担。

Swift 输入依次为 round032/exceptional_signed_mask_result.json、round033/mm_completion_result.json、round034/macbook_s32_flat_audit_result.json 与输出路径。可用 swiftc -O 编译 R35LocalDimensionMetal.swift，二进制与模块缓存使用该机授权的临时目录。MacBook 曾因 SSH 横幅超时延迟，最终运行和编译均退出0。
