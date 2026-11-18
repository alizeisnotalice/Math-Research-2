# 第53轮复现

运行环境为 ARM64 macOS、Python 3.15.0、Apple clang 21.0.0 和现有 Swift Metal 后端。
以下命令从本证书目录运行。用 data 卷的空目录保存输出，避免覆盖归档。

## 解析与真条件词的整数核验

```sh
python3 -B analysis/audit_witness_prefix.py | cmp - analysis/audit_witness_prefix.stdout
python3 -B geometry/audit_distance_layers.py --workers 12 | cmp - geometry/distance-layers.json
python3 -B geometry/verify_contact.py \
  ../../computation/results/exact-r001.json \
  ../../computation/results/d4-generic-n5-r005.json | cmp - geometry/contact-verification.json
```

重算首接触表时，将 geometry/contact_geometry.cpp 和 run_contact.py
复制至新目录。用共享 clang 缓存，以 `clang++ -O3 -std=c++20 -Wall -Wextra -Werror`
编译为 contact_geometry。三、四阶分别用命令 `./contact_geometry 3 3` 及
`./contact_geometry 4 3`；五阶用 `python3 -B run_contact.py n5-edge 0 3`
和 `python3 -B run_contact.py n5-local 3 6`。五阶前半在 MacBook 复算后逐字节
比对，本档不另存重复 TSV，远端回执的哈希可用 n5-edge.tsv 校验。
编译产物哈希、测时、内存和 PID 可以变化。

## 联合缝峰的十二阶分母

在新目录复制 seam 的 joint_seam.cpp、run_joint.py、merge_joint.py、
pratt-61bit.json 和 n12-plan.json。按以上 clang 参数编译 joint_seam。
运行 `python3 -B run_joint.py batch`，再运行
`python3 -B merge_joint.py batch upper.json --prior <课题目录>/certificates/r052/t5`。
除了工人回执的哈希，独立重编译的科学账本与 joint-bound-final.json 完全一致。
原始构造分子从第52轮归档读取。精确联合过滤的解释见
[联合证书](seam/JOINT_SEAM_UPPER.md)。

## 本机 Metal

两个 Swift 源使用同一份冻结的 GPU 内核，只增加宿主端非随机归约。
contact.swift 可由 `python3 -B make_contact.py transport.swift contact.swift`
生成。用 `swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift`
分别编译为 transport 与 contact；复制 run_metal.py 至同一目录。
实际依次运行 gate、replay、n24、contact-gate、contact-n16 五个模式。
末两个模式在 n24 退出后开始。参数、种子域、运行终态与数据哈希见各模式回执。

```sh
python3 -B metal/verify_metal.py --r052 ../r052/local-metal/control/n16.jsonl \
  | cmp - metal/verification.json
```

三阶17粒子跟踪使用独立 Python 图程序逐粒子重算两种宿主归约。
十六阶两次复算与第52轮同配置下的所有旧第一合并矩、峰距直方图和
对数归一化逐重复吻合。首接触可达性使用整条前缀和准确的插槽权 \(L=j-1\)；
只依两个峰位置的角点路径成本使用第一合并权。两者是不同的数值日程。
本轮采样只用于有限探路，没有混合误差界。
