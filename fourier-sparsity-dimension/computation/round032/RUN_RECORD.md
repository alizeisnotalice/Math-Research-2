# 第32轮实际运行记录

主 Agent 在最终归档目录运行所有脚本，均exit 0。Python 3.15.0、macOS 27.2 arm64；C为Apple clang 21.0.0，-O3 -std=c11 -Wall -Wextra，编译无警告。全部是GF(2)、整数FWHT或位计数，无随机种子、无浮点误差容限。时间字段只用于调度，不参与数学判断。

- `rm48_weight32_kernel_audit.json`：0.547405秒（C字段为进程CPU时间，其余为墙钟）。
- `transvection_avoidance_result.json`：1.929230秒（C字段为进程CPU时间，其余为墙钟）。
- `k5-transvection-candidates.json`：1.284496秒（C字段为进程CPU时间，其余为墙钟）。
- `r18-tensor10-transvection-histogram.json`：1.050912秒（C字段为进程CPU时间，其余为墙钟）。
- `exceptional_signed_mask_result.json`：1.804388秒（C字段为进程CPU时间，其余为墙钟）。
- `signed_plateaued_enum_result.json`：0.006327秒（C字段为进程CPU时间，其余为墙钟）。

9类代表各由独立进程完整枚举97155个RREF三维子空间；最终输入移除机器专属clone路径，脚本改为相对路径后重跑。另有主Agent集合实现与agent位集实现逐项一致的前置交叉核验。完整分类本身来自已发表来源，本轮未重跑作者分类程序。

几何脚本检查所有12870个F₂⁴八点集的避差核降维等价，另验证一个实际八变量bent剪切实例。早期临时版本曾错误地用加法合并重复差分位，已改为按位OR；归档只保留修正后实际重跑并通过的源码。

两个十变量种子各查522753个非平凡可逆剪切；均无重量64差分。所有导数实际检查平衡，种子实际检查bent；搜索范围不能外推。第18轮seed从相对路径输入读取，输入JSON记录其哈希。

唯一外包负类的符号由两种枚举独立穷尽1470943个规范状态，均恰1解。Python与C的支撑顺序不同，验证器按负号点集合比较，且重新直接求128个Walsh系数验证16个非零值均为±16。C结果的source_sha256与compiler字段由主Agent在运行后补入；复跑原C stdout没有这两项，可按源码哈希补入，不把这些元数据当成数学证书。

本轮本机9路精确核验批次很短，完成后均退出；用户要求查看满载时本机连续实测0% idle，但包含其他课题负载，不归因于A28独占。本轮未使用MacBook、GPU或MATLAB，不以历史运行冒充本轮计算。大内存任务没有启动。

复现（仓库根目录，C二进制置于data临时区）：

```sh
python3 fourier-sparsity-dimension/computation/round032/audit_rm48_kernels.py
python3 fourier-sparsity-dimension/computation/round032/check_transvection_avoidance.py
python3 fourier-sparsity-dimension/computation/round032/k5_transvection_candidate_count.py
python3 fourier-sparsity-dimension/computation/round032/r18_tensor10_transvection_histogram.py
python3 fourier-sparsity-dimension/computation/round032/check_exceptional_signed_mask.py
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang -O3 -std=c11 -Wall -Wextra fourier-sparsity-dimension/computation/round032/signed_plateaued_enum.c -o /Volumes/data/CodexTemp/T2/a28/r032-signed/signed_plateaued_enum
/Volumes/data/CodexTemp/T2/a28/r032-signed/signed_plateaued_enum
python3 fourier-sparsity-dimension/computation/round032/verify_artifacts.py
```

验证器核对输出及小证书，没有重跑完整外部分类或重新完成所有枚举。原始数学条件的证明见proofs/032-transvection-rm-obstruction-zh.md，manifest仅检查归档完整性。

继承材料（没有在本轮重新执行的历史证明/证书）：

- `proofs/019b-six-variable-projection-audit-zh.md`：`767963dbc5de083601776dc2a3fb957fb25fb56efa725bfbc0b78a2fd54e42f3`。
- `proofs/019-projectable-double-blocks-zh.md`：`b66a34ddd5cda30c8cb7aa2ca57030b6b44b0f80a202aee773d34e8594650495`。
- `computation/round018/graph-cases.json`：`52d992aa77f2075575fee144069602babc54266e03b7ee3deec48ed9ba9daa81`。
