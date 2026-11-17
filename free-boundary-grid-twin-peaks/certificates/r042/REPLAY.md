# 弯曲切割精确证书复现

源码自写，C++20/Apple clang21，Python3.15标准库。n5中央二值子族生产152个轨道，
分别代表32个对角事件和496个无序不同事件对。n5-ledger.json保留每个代表计数、
完整轨道成员及原始输出SHA。原始raw/时间/逐child收据完整本地封存；公开提交
合并后的最小整数账本，不批量提交152份资源日志。evidence保留最终退出与聚合回执。
此子族使用原无剪枝版本。扩大178族另用已核验的严格零剪枝版本，4134个轨道全部实际完成。扩大账本n5-lipschitz-ledger.json省去冗余members；审计从178曲线重建全部15931个无序对，逐代表核验两个乘数。

将PKG指向本目录，OUT为全新仓外目录：

```sh
mkdir -p "$OUT" /Volumes/data/CodexTemp/cache/clang
env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
 clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pedantic \
 "$PKG/source/curved_cut_dp.cpp" -o "$OUT/curved_cut_dp"
python3 "$PKG/source/audit_small.py" "$OUT/curved_cut_dp" "$OUT/n3-exhaustive"
python3 "$PKG/source/audit_certificate.py"
python3 "$PKG/source/audit_diagonal.py"
python3 "$PKG/source/audit_diagonal.py" --lipschitz
python3 "$PKG/source/audit_lipschitz_certificate.py"
python3 "$PKG/source/make_jobs.py" 5 "$OUT/jobs-n5.json"
python3 "$PKG/source/make_lipschitz_jobs.py" "$OUT/jobs-n5-lipschitz.json"
```

n3独立枚举全部9!标号及36交集。改进版另外复查横竖16事件的136个无序交集：

```sh
env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
 clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pedantic \
 "$PKG/source/curved_cut_dp_pruned.cpp" -o "$OUT/curved_cut_dp_pruned"
python3 "$PKG/source/audit_small_pruned.py" "$OUT/curved_cut_dp_pruned" "$OUT/n3-pruned"
```

n5完整重算逐个使用相应生成清单中的代表，原版例：

```sh
"$OUT/curved_cut_dp" 5 "$OUT/n5-22222.json" 22222 22222
"$OUT/curved_cut_dp" 5 "$OUT/n5-22222-33333.json" 22222 33333
```

一worker核心512MiB，全部152代表需依次或按真实可用内存调度；历史生产最多8worker，
单job最大RSS538,525,696字节。完整生产已经完成、wait/reap，controller/caffeinate及
任务进程已消失。小n父Agent重放不冒充第二次完整n5生产。
178族4134代表采用curved_cut_dp_pruned：1051个代表因严格几何必要条件返回0且不分配大表，其他执行完整递推；最大RSS538,542,080字节，也已完成并核验本任务PID消失。

原直切程序额外交叉n4=41328967680，n5两个直切各301412212232374763520、
其交集14001081747604852400；曲线DP逐项一致。
此外audit_diagonal独立对每侧诱导图做整数子集递推，再乘二项交错数，
分别逐个核对32/178个n5单事件计数及其和；它不使用全图多约束递推。
