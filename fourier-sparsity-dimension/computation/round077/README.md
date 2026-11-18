# 第77轮复现

partition.cpp 使用 C++20，以包含点 0 的首块位掩码除二取整作为分片号，递归时各块按最小元素排序。编译命令为 clang++ -O3 -std=c++20 partition.cpp -o partition。两台设备分别运行 ./partition 0 8192 15 及 ./partition 8192 32768 12。四项判定见 test 函数。分片总数等于 Stirling 第二类数 S(16,1) 至 S(16,4) 之和。编译环境为 Apple clang 21，在 Apple Silicon 上运行。

在本目录用 python3 -B 分别运行 gr-verify.py、verify_tree.py、verify_aff.py、verify_odd.py、check_flag_rounding.py、check_odd_walsh_hierarchy.py。运行 check_odd_walsh_mod8.py 时指定 exhaustive 或 sample。运行 check_odd_walsh_metal.py 需要 macOS Metal。整数及有理数判定的主要输出和报告包含在同目录。抽样的种子、环境、维数及检查数保存在各 JSON 中。CPU 线程数为启动配置，本轮没有测量划分作业的持续核心利用率。
