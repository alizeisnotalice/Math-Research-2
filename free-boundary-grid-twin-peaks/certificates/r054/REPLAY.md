# 第54轮复现

本轮在 ARM64 macOS 使用 Python 3.15.0、Apple clang 21.0.0 与已有 Swift Metal 程序。
以下从本证书目录执行，SCRATCH 为 data 卷中的空目录。

## 早接触与标记路径的整数核验

```sh
python3 -B transport/verify_markers.py ../r053/geometry | cmp - transport/verification.json
python3 -B transport/verify_trapped.py ../r053/geometry | cmp - transport/trapped.json
python3 -B bridge/verify.py | cmp - <(sed '1d;$d' bridge/verify.txt)
```

重算 C++ 路径表时，用 `make_census.py ../r053/geometry/contact_geometry.cpp`
加输出路径生成 marked_paths.cpp，或直接使用归档的同一源码。
使用 `clang++ -std=c++20 -O3 -Wall -Wextra -Werror` 及共享 clang 模块缓存编译。
第一个参数为边长、第二个为线程上限，归档表的边长为3、4、5。
计算限于首接触秩不超过7及长度至多4的简单路径。
近区按格点奇偶分别作标记，结合两张输出容量时各取一半权。
五阶在 MacBook 重新编译并运行六个轨道，结果逐字节一致，另存运行和哈希回执。

## 固定词界面交错

```sh
python3 -B interface/interface_weights.py "$SCRATCH/interface.json"
cmp "$SCRATCH/interface.json" interface/interface-gate.json
```

对小矩形完整织合按真实峰集校验，对二色反射及多色序锥作整数枚举，
并用较大格中构造的固定局部词计算完整交错权。

## 24阶 Metal

将第53轮的 contact.swift 与本轮 run_metal.py 复制至空工作目录，
用 `swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift` 编译为 contact。
运行 `python3 -B run_metal.py contact-n24`，得到原始运行目录。
本次实际复用了同哈希本机二进制，只改变运行配置，种子域为9e54100000000001。

```sh
python3 -B metal/verify_n24.py metal ../r053/metal/n24/run.jsonl \
  ../r053/metal/contact.swift | cmp - metal/verification.json
```

检查逐重复早期秩和原二维距离整数边缘，使用对数归一化聚合。
第53轮比较使用另一种子域和更多移动，程序直接从归档读取。
两份均为有限诊断，未证明混合精度。
