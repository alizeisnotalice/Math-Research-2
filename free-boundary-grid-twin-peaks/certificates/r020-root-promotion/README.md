# 第20轮：最高根与两最高根的精确计数

三个Markdown文件给出解析证明。`root_deletion_identity_n3.py`、`two_root_block_identity_n3.py`分别完整枚举9!标号核验两恒等式。`two_root_block_ie.cpp`包含仓库相对路径的既有独立集DP，交叉n3并计算n4；`verify_two_root_block.py`读取结果核对。主Agent已重编译并重跑n4，输出JSON与归档完全相同。

在本目录启动，将新产物写到临时盘：

```sh
TASK_TMP=/Volumes/data/CodexTemp/T2/a022-r020-promotion-replay
mkdir -p "$TASK_TMP"
c++ -O3 -std=c++20 -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang-modules two_root_block_ie.cpp -o "$TASK_TMP/two_root"
"$TASK_TMP/two_root" 4 > "$TASK_TMP/n4.json"
python3 root_deletion_identity_n3.py > "$TASK_TMP/root-n3.json"
python3 two_root_block_identity_n3.py > "$TASK_TMP/two-root-n3.json"
python3 verify_two_root_block.py
```

继承DP的main函数被宏改名，编译器会提示该未调用函数缺末尾return；实际新main独立返回，核验不调用该入口。源及完整账本记录于本轮manifest。

最强新结论只涉及局部2至4峰标号的总质量，不能把它替换为双峰计数。删峰的平均逆像重数仍未控制。
