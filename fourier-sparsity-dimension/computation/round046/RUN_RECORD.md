# 第46轮计算记录

本轮计算使用精确有限域消元、整数 Walsh 变换和有理数概率。通用树定理由证明文件中的逐路径消元、频谱隔离和仿射陪集体积承担。

|程序及结果|实际运行和核验范围|
|---|---|
|[混合树主程序](mixed_parity_audit.py) · [结果](mixed_parity_result.json)|本机 Python 3.15.0，6 个进程，七个预设模板，给定叶赋值集的完整枚举，无随机种子。178088 组叶赋值，2198428 个根 Walsh 系数，1232544 个隔离后系数。归档检查复跑一次，除计时外所有输出字段一致。|
|[商空间改写](verify_mixed.py) · [结果](verify_mixed.json)|本机 Python 3.15.0，确定性种子 0 至 399。选择空间 3 个私位、2 个共享位、1 个独立数据位。原树与改写树在每个输入给出相同真值；对 25600 个联合频率位置及 3062 个标记私节点进行检查。1353 个标记节点用到共享分叉。复跑输出逐字节一致。|
|[共享偏移核验](audit_mixed_offset_isolation.py) · [结果](audit_mixed_offset_isolation.json) · [独立反例检查](nonaffine_private_twist.py) · [结果](nonaffine_private_twist_result.json)|本机 Python 3.15.0；第一份代码覆盖共享维数 2 至 6 的五个反例及五个仿射剪切例，第二份直接确认维数 2、3 的完整谱，并核验一个不付非线性代价的全局界失败。两份均无随机步骤，归档时独立复跑输出逐字节一致。|
|[门秩取等程序](joint_gate_equality.py) · [结果](joint_gate_equality.json)|本机 Python 3.15.0，$m=0,\ldots,8$；逐地址和全输入计算 Walsh 谱并以 Fraction 核对门质量、加权下界及形式节点数。|
|[Swift 与 Metal 源码](R46MixedParityBlockerMetal.swift) · [本机结果](local_r46_mixed_parity_result.json) · [MacBook 结果](macbook_r46_mixed_parity_result.json) · [MacBook 运行回执](macbook_r46_mixed_parity_receipt.txt)|本机 Apple M6 和 MacBook Apple M5 Pro，Apple Swift 6.4.2。每机 12 棵指定树和 192 棵固定种子树（十六进制种子 a280460edabc2345），33656 个整数 Walsh 系数逐项匹配，22612 次屏蔽后系数检查，3606 次门容量检查，943 个标记节点。两结果仅设备名称及 GPU 计时不同。MacBook 回执所列源码与结果文件 SHA-256 已核对；本机独立编译运行后确认该核验程序已退出。|

MacBook 运行脚本 [run_r46.sh](run_r46.sh) 是当次远端工作区使用的原始脚本。归档文件中不含可执行程序、缓存或编译中间件。MacBook 当次回执所列 Swift 源码 SHA-256 为 6f24664314e9d03f32de8048fa4f71d000956f819445eb822c689f40960514d8，结果为 7de58596b3651305a31c4130f627730aefd295c8007730805c0bb9672923eec6。

从仓库的本轮计算目录运行下列命令，临时输出放在 data 卷。

```sh
out=/Volumes/data/CodexTemp/T2/a28/r046-replay
mkdir -p "$out"
helper=../round045/parity_tree_audit.py
python3 -B mixed_parity_audit.py "$helper" "$out/mixed.json" --workers 6
python3 -B verify_mixed.py > "$out/quotient.json"
python3 -B audit_mixed_offset_isolation.py > "$out/offset.json"
python3 -B nonaffine_private_twist.py "$helper" "$out/twist.json"
python3 -B joint_gate_equality.py > "$out/equality.json"
mkdir -p /Volumes/data/CodexTemp/cache/swift-module-cache /Volumes/data/CodexTemp/cache/clang-module-cache
export TMPDIR="$out" TMP="$out" TEMP="$out"
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang-module-cache
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-module-cache R46MixedParityBlockerMetal.swift -o "$out/r46_metal"
"$out/r46_metal" "$out/metal.json"
```

主枚举及独立非仿射反例程序复用第45轮的 GF(2) 和 Walsh 辅助函数。本轮实用辅助文件 SHA-256 为 6dfb101601149eb93039d31a494d08cd02cfb6f3f37fad6b9c5ae86d9ae6a9a2。主枚举结果包含耗时；Metal 计时及设备名称依运行环境而异。
