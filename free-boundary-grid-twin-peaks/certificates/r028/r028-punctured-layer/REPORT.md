# r028：选定最高点孔洞的模数分层容斥

状态：有限精确计算的实现与小规模交叉核验已完成；n=7 生产运行不属于本目录的本次小规模核验，且不得据此作任何渐近结论。

## 计数对象和递推

在 `P_n square P_n` 固定全局最高标号位置 `v`。记 `D=N[v]`，容许的指定峰集
`I` 是 `V \ D` 的独立集。以原图闭邻域记 `C(I)=N[I]`，则固定 `v` 后的补全计数为

```
C_v(empty) = (N-1)!
C_v(I) = sum_{u in I} C_v(I\{u}) / |N[I]|.
```

因为 `I` 避开 `N[v]`，`v` 不在 `N[I]`，每个分母均在 `1,...,N-1`。程序逐层枚举这些独立集并以四个模数分别计算：

```
A2_v = sum_{I != empty} (-1)^(|I|-1) |I| C_v(I)
B_v = 1/(N-1) sum_I (-1)^|I| C_v(I) |V \ (N[v] union N[I])| C_v(I).
```

这是有限组合恒等式的程序实现；公式本身另由主线核验。本实现不把容许峰集的闭邻域错误地要求避开 `N[v]`：只断言它不含顶点 `v`。

## 模数、恢复和资源安全

四个经 r027 Pratt 认证的素数为

```
9223372036854775783, 9223372036854775643,
9223372036854775549, 9223372036854775507.
```

其乘积为 252 bits，而 n=7 时每个非负补全计数至多 `(49-1)!`，仅 203 bits。因此四个残数的 CRT 恢复唯一。每个素数都大于 48，故递推分母在模域可逆。

值表按大小只保留相邻两层；前层的 64-bit mask 到索引哈希表采用容量为最少的二次幂且装载率不超过 5/8。中心孔 `n=7,v=24` 的精确层统计为 293,815,906 个独立集，最大转移是 `k=12`：前层 56,669,220、当前层 52,506,652、哈希槽 `2^27`。模型 payload 为

```
8*(56,669,220+52,506,652) + 12*2^27 = 2,484,019,712 bytes = 2.3134 GiB.
```

这是受 7 GiB 程序内 payload guard 约束的模型数，不是 RSS 上界，更不能替代运行前整机内存检查。n=7 生产只可四个模数顺序运行。

## 已实际运行的小规模核验

冻结 C++ 源码：`punctured_layer_mod.cpp` SHA-256
`5f93fa39695e0f4badabcbabaad2befab920fa1590ea2337fe47a09201587db9`；二进制 SHA-256
`e57308e1d422b949ed43c39f8939f3fb046d7041e6f1adf61a93d48b7994c37c`。
编译器为 Apple clang 21.0.0，Python 3.15.0。64 位位掩码边界测试通过（含 bit 63）。

单线程逐个运行 72 个 checkpoint：`n=3,4,5,6` 的每个 D4 代表位置、每位置四素数。每个 checkpoint 具有独立的 stdout、`/usr/bin/time -l` 资源记录和 `exit_status=0`。`verify_residues.py` 以 Python 整数 CRT 重构后，18 个 `(n,v)` 均匹配相应精确表：

| n | 核验的 v | 表来源 | 结果 |
| --- | --- | --- | --- |
| 3 | 0, 1, 4 | r027 punctured-table-n3 | 3/3 精确匹配 |
| 4 | 0, 1, 5 | r027 punctured-table-n4；另有独立 subset-birth oracle | 3/3 精确匹配 |
| 5 | 0, 1, 2, 6, 7, 12 | r027 punctured-table-n5 | 6/6 精确匹配 |
| 6 | 0, 1, 2, 7, 8, 14 | r028 punctured-table-n6-assert | 6/6 精确匹配 |

一个实际最重的小检查 `n=6,v=14,p=0,T=1` 用 0.26 s wall、19,743,128 B peak memory footprint、0 swaps；它不预测 n=7 的 RSS。完整整数结果和来源哈希在 `finite-crosscheck.json`；校验命令及自身哈希在 `SHA256SUMS`。

## 复现

```
env CLANG_MODULE_CACHE_PATH="$CLANG_MODULE_CACHE_PATH" \
  clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pthread punctured_layer_mod.cpp -o /path/to/punctured_layer_mod
/path/to/punctured_layer_mod --boundary-test /path/to/boundary-test.json
/path/to/punctured_layer_mod 6 14 0 1 /path/to/n6-v14-p0.json
python3 verify_residues.py --runs runs \
  --table-root "$A022_R027_CERT_DIR/r027-global-normalizer" \
  --table-n6 "$A022_R028_CERT_DIR/r028-punctured-n6/punctured-table-n6-assert.json" \
  --output finite-crosscheck.json
```

本轮未自行启动 n=7，也未使用 GPU：这是精确模整数/哈希计算，CPU 路径保持整数算术。小核验进程已退出；生产进程由主线独立管理。
