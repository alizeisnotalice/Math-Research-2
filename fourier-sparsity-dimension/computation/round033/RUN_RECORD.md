# 第33轮实际运行记录

## 本机发现与精确验证

最终源在归档目录实际运行，Python 3.15.0、macOS27.2 arm64，全部exit0，无随机种子。匹配搜索用精确GF(2)与16×16二分图，在第28个坐标拆分得到完美匹配，约0.005390秒；不声称完成70种拆分或全部MM/bent分类。物理及频率约束来自第32轮唯一规范signed mask。

独立脚本从两个256位真值重建函数，逐点检查MM表、置换性与h公式，重算全部Walsh、ANF和物理差集，约0.001375秒。两函数全谱均±16，32点差块的差集恰为128点超平面；对偶32点块的182点差集避开span{97,98,20}的全部七个非零向量。

张量/剪切脚本约0.192489秒，检查r=0..3：张量维8/10/12/14，另两变量线性剪切维10/12/14/16。全部构造及其翻号均bent，重量和剪切差分公式逐点通过。一般r结论由解析Walsh求和及维数交公式证明，不由有限检查外推。

## MacBook CPU 与 GPU

MacBook通过已配置SSH实际执行。Apple M5 Pro、Metal4、Swift6.4.2；最终Swift源SHA256 `cce151ef81d0fa81e24793139f4b7eaa410c37ba6965653c51acc4c07549a67e`。CPU为Int32迭代FWHT，GPU为Metal Int32全局分阶段双缓冲FWHT，无浮点舍入，最大域65536点，所有中间绝对值≤65536，远小于Int32范围。

实际运行只核验张量对 b+u·v 与 b+u·v+h1_{v=0}，r=0..4、维8/10/12/14/16，共五组十个函数。所有系数CPU/GPU逐项相同，bent违反数均零，差分重量32/64/128/256/512。**远端没有运行另加t,s的线性剪切式**；后者由本机上述脚本核验。避免把张量翻号称为已运行的线性剪切实验。

编译与运行exit0，含编译1.434471秒，程序内0.010222秒。早期一次Swift空格编译错误、一次未允许hex前导零省略的解析断言失败，已修正后完整重跑；失败运行不作证据。进程已退出，未启动随机补全搜索、MILP长任务或MATLAB。此次短任务不声称持续满载。

远端输入是发现时result.json，归档输入则为最终相对路径脚本重跑结果，两份solution对象已逐字段核对相同。归档简化回执记录原始输入哈希及共同solution规范哈希；没有把不同元数据/计时的文件说成逐字相同。复现时直接使用归档mm_completion_result.json，输出输入哈希将相应改变。

## 复现

在仓库根目录：

```sh
python3 fourier-sparsity-dimension/computation/round033/search_coordinate_mm.py
python3 fourier-sparsity-dimension/computation/round033/independent_mm_witness_audit.py
python3 fourier-sparsity-dimension/computation/round033/check_tensor_shear.py
python3 fourier-sparsity-dimension/computation/round033/verify_artifacts.py
```

Metal源接受三个文件参数，可复制到任何已授权且Metal可用的Mac，用该机已确认的临时路径编译。以本机data缓存为例（此为复现命令，实际GPU运行在MacBook）：

```sh
mkdir -p /Volumes/data/CodexTemp/T2/a28/r033 /Volumes/data/CodexTemp/cache/swift-module-cache
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-module-cache fourier-sparsity-dimension/computation/round033/R33MMTensorMetal.swift -o /Volumes/data/CodexTemp/T2/a28/r033/check_metal
/Volumes/data/CodexTemp/T2/a28/r033/check_metal fourier-sparsity-dimension/computation/round033/mm_completion_result.json fourier-sparsity-dimension/computation/round032/exceptional_signed_mask_result.json /Volumes/data/CodexTemp/T2/a28/r033/reproduced_metal.json
```

验证器重新执行种子整数FWHT和小核证书，并核对运行产物及绑定；不会远程重启GPU，也不代替一般解析证明。manifest仅用于文件完整性。
