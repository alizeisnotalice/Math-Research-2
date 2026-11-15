# final-r005 来源、版本与复现范围

本目录是第 5 轮 beta-support 精确计数的紧凑归档；不包含编译二进制或逐轨 stderr。所有归档文件的 SHA-256 在 `SHA256SUMS`。

认证执行使用 `beta_zero_layered_orbit_certified.cpp`，SHA-256 为
`e0c769ad793d63391dc3e853b7170680afa04dfde4b53c1c656f9d1fd875c3b5`。
它由 Apple clang 21.0.0 以 C++17、`-O3 -pthread` 编译。运行前/后唯一命名二进制的 SHA-256 都是
`945b1cb671f8c04bc36c5e267de79f1e3c73783644f4b2e3ef5a7ae88d622cb3`；二进制按归档要求不复制。

`run_n7_beta_orbits_certified.py` 使 12 个 D4 无向边轨道依次以独立进程执行，每进程 12 workers，且仅在前一轨进程退出、JSON 解析通过后启动下一轨。`orbits/manifest.json` 保存准确命令、每轨 UTC 起止时刻、exit code、运行前/后源码和二进制哈希；`orbits/orbit-*.json` 和 `orbits/orbit-*.time` 分别是逐轨计数输出和 `/usr/bin/time -l` 资源记录。

`verify_beta_zero_layered.py` 仅使用 Python 3.15.0 标准库。它核验 n=3/5/6 基准、历史完整 n=7 JSON，及分轨回放的 12 个 exit code、资源文件、覆盖数和加权总和。实际输出证书为 `beta-zero-layered-orbit-certified-verification.json`。

此归档的 n=7 结论是有限精确恒等式
\[
B_0(7)=131529421404957780122938007473728844832061946591752192.
\]
它不估计 \(\mathbb E_2\beta\)，不推出 \(R_n\) 或原题的渐近结论。

历史的单进程完整 n=7 运行二进制哈希未在覆盖前保存；本目录不倒填或声称恢复该哈希。这里的分轨认证是独立、带源码—二进制哈希对应的新运行，且数值与历史完整 JSON 完全相等。
