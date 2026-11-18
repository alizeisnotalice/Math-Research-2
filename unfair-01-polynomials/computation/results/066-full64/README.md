# 64 次第 7 项缺口的完整精确排除

[完整归约与证明](../../../proofs/066-c7-forced-bits-to-30-and-degree64.md)。

第 0 至 35 位除了第 7 位均为 1。第 38 位和首项为 1，第 60 至 63 位为 0。其余 23 位任意。两套实现分别从模 163 的普通余式和模 1279 的负倒多项式伪余式计算必要的导数 gcd 门槛。两个首用素数的 8,388,608 个目标均无通过项。独立核验器检查全部 256 个分片，用 SymPy 复算 68 个目标在两素数下的精确 gcd 次数。两个实现均以 \((z^{32}-1)^2\) 作正命中回归。

使用 Apple clang 21.0.0、Python 3.13.18、SymPy 1.14.0。在本目录重放如下，输出写到包外尚未存在的目录。

    root=/Volumes/data/CodexTemp/T2/a30-full64-replay
    mkdir -p "$root/bin"
    clang++ -std=c++20 -O3 scan64.cpp -o "$root/bin/scan64"
    clang++ -std=c++20 -O3 scan64_reciprocal.cpp -o "$root/bin/scan64_reciprocal"
    python -B dispatch64.py "$root/bin/scan64" "$root/ordinary" --workers 10
    python -B dispatch64.py "$root/bin/scan64_reciprocal" "$root/reciprocal" --workers 10
    python -B verify64.py "$root/ordinary" "$root/reciprocal" "$root/bin" "$root/fixtures" "$root/verification.json"

[平方次数七边界](boundary/README.md)由不同的参数遍历独立闭合。MacBook 完成旧的 162 参数搜索，每半枚举 \(2^{26}\) 个码，各参数无命中，回执为 mac-boundary-verification.json。普通和负倒 gcd 的完整结果同时覆盖较大平方次数和有理重因子。
