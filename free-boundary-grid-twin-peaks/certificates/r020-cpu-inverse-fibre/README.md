# 第20轮：一峰逆纤维 CPU 精确测量与 SMC

`inverse_fibre_smc.cpp` 的源SHA-256为 `e6eff498c7274007e2417414f87eafcaf1e5b3ecf8d02edd8c0a835fd6be2289`；差异文件给出相对r019根矩核的全部修改。随机数/整数CDF/MH保持原核，峰数约束改为onepeak，终端m2/m3使用整数DP；每个终端词另按rank与邻点检查唯一峰。此程序直接估计A1计数，和GPU的概率Z相差N!。

全部29,568个3×3 onepeak词的C++ m2/m3表，与独立Python真实提升/DP逐字节相同，sum m2=63,552、sum m3=64,512。只存一份CSV，另一个可由`make_n3_fibre_table.py`重建；相等回执见 `exact-table-comparison.json`。主Agent另行重编译并运行全表审计，见 `root-table-replay.json`。

n3 R512的10%门槛为探索后选用的数值诊断，不是预注册测试，不认证概率误差。n7用精确计数校准；n10、n11是有限数值证据。主Agent `summarize_inverse.py` 用Decimal80检查完整rep覆盖、种子、提议次数和加权比，再给跨rep权重ESS/最大权重与delta-SE。它不把particle parent ESS当作跨rep ESS。

n11两批用不同seed_base，实际为8线程R256与4线程R128，并行使用12线程；两批已完成、退出，原始CSV/metadata/资源分别保存。两个domain单独报告，未伪装成同seed增加独立样本。所有数字不承担渐近证明或混合认证。

复现时在本目录启动，编译/输出放临时盘：

```sh
TASK_TMP=/Volumes/data/CodexTemp/T2/a022-cpu-inverse-replay
mkdir -p "$TASK_TMP"
c++ -O3 -std=c++20 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang-modules inverse_fibre_smc.cpp -o "$TASK_TMP/inverse"
"$TASK_TMP/inverse" audit-n3 "$TASK_TMP/n3-cpp-fibres.csv"
"$TASK_TMP/inverse" 11 65536 256 16 8 "$TASK_TMP/n11.csv" 0xA0222B00
python3 summarize_inverse.py n11-p65536-r256-l16.csv n11-p65536-r128-l16.csv --output "$TASK_TMP/n11-summary.json"
```

`verify_n3_calibration.py` 的完整双表检查需要先在临时副本中运行 `make_n3_fibre_table.py` 生成另一张表；已归档相等证据不替代重新生成。该辅助脚本引用相邻 `r020-inverse-fibre/exact_inverse_fibre_n3.py`。
