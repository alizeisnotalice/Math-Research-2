# A28：签名选择器族的完整精确枚举

对象为
\[
g(u,y)=\sigma(u)(-1)^{b(u)\cdot y},\qquad
u\in\mathbb F_2^3,\quad y\in\mathbb F_2^2,
\]
其中 \(b(u)\in\mathbb F_2^2\)，\(\sigma(u)\in\{-1,1\}\)。整体乘以 \(-1\) 不改变 Fourier 支撑，故枚举规范化为 \(\sigma(0)=1\)：共有 \(4\cdot8^7=8{,}388{,}608\) 个实例。对每个 \((a,c)\)，程序以整数计算
\[
8\widehat g(a,c)=\sum_{u:b(u)=c}\sigma(u)(-1)^{a\cdot u},
\]
再以 \(\mathbb F_2\) 消元求支撑线性秩。没有浮点零阈值。

MacBook 的实际运行使用 `clang++ -O3 -std=c++17`，单低优先级 CPU 进程完整扫描耗时 `0.921589` 秒；原运行结果位于临时目录，最终相同内容已归档到本目录 remote-full.json：
`/Volumes/data/CodexTemp/T2/a28/remote-selector/remote-full.json`。
输出 SHA-256 为
`a9b11410d60c4f17f339cb50284d0b020e1733b153fd9e4a29fd5b9a2b1bfdcb`。

完整范围内，每个出现的谱秩的最小稀疏度为：

| d | 最小 s | s/d² |
| --- | ---: | ---: |
| 1 | 1 | 1 |
| 2 | 4 | 1 |
| 3 | 4 | 4/9 |
| 4 | 8 | 1/2 |
| 5 | 10 | 2/5 |

`verify_selector_result.py` 在本机重新以整数 Walsh 变换和独立的 \(\mathbb F_2\) 消元验证各维返回极小值候选；收据位于
`/Volumes/data/CodexTemp/T2/a28/remote-selector/local-reverify-remote-full.json`。
这是一族有限构造的穷尽计算，不能证明一般二次下界，也不能构成无穷反例族。

本机 GPU 实测使用 PyTorch 2.15.0 / MPS，批量 65,536、三次重复；CPU 和 MPS 的计数相同，MPS 中位时间 `0.01868025` 秒，CPU 中位时间约 `0.006` 秒。因此精确扫描使用 CPU；校准记录在
`/Volumes/data/CodexTemp/T2/a28/remote-selector/local-cpu-mps-calibration.json`。

复现命令：

```sh
clang++ -O3 -std=c++17 selector_search.cpp -o selector_search
./selector_search 0 8388608 remote-full.json
python3 verify_selector_result.py remote-full.json --out receipt.json
```

主Agent整合补充：execution-summary.json记录硬件、编译器、源码哈希和本机完整重跑；parent-verification.json通过完整真值表字符求和、ID及覆盖计数检查；parent-cpu-mps-calibration.json逐项比较全部Walsh输出。旧校准只比较合计数，不应描述为逐项验证。新版脚本的结果与旧版校准记录分别保存。复现编译时使用/Volumes/data/CodexTemp/cache/clang模块缓存，可执行文件放CodexTemp。

远端工作目录：`/Users/zhengzhihao/CodexRemoteMath/a28-selector`。实际命令：`clang++ -O3 -std=c++17 selector_search.cpp -o selector_search`，随后 `nice -n 15 ./selector_search 0 8388608 remote-full.json`；验证用 `/usr/bin/python3 verify_selector_result.py remote-full.json --out remote-full-receipt.json`。远端源码SHA-256与本地相同，见execution-summary.json。Apple clang 21.0.0为主Agent在运行结束后探测所得，原编译输出没有版本行。
