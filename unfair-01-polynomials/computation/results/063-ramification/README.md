# 第 7 项缺口的分歧 Hankel 证书

[任意次数证明](../../../proofs/063-c7-ramified-hankel.md)。对 163、1279 两素数，低端幂和的七项正规化差商分别有非零的四阶 Hankel 行列式 111、217。模平方自由部分与导数的 gcd 次数至多 3 时，惯性群的偶数阶只可能是 2、4、6。阶 2、4 给不足的 Hankel 秩；阶 6 给一个移动六次根的单换位，与首三项幂和矛盾。

Python 3.13.18 与 SymPy 1.14.0。重现行列式和有限簇类型，并逐项复核 12 个具体目标的两个 gcd 都是一次式。

    python3 -B verify_ramified_hankel.py --out /Volumes/data/CodexTemp/T2/a30-ramified-hankel-new.json
    python3 -B verify_ramified_collision.py inputs-gap7-degree56.json --out /Volumes/data/CodexTemp/T2/a30-ramified-collision-new.json
    python3 -B verify_ramified_packet.py --out /Volumes/data/CodexTemp/T2/a30-ramified-packet-new.json

C++20 中间相遇生成 425 个字。生成时在模 163 的 \(-2\) 和模 1279 的 \(2\) 同时强制值与导数为零，再由 Python 用有理域平方自由检查、剩余域次数及数值根准备筛到 12 项。未处理的字没有排除状态。前一道筛选未覆盖的重数型现在由 Hankel 定理排除。六个已运行 Metal 分段保存并退出，后续不再启动。
