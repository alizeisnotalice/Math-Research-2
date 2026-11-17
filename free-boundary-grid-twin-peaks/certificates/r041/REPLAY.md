# 第41轮复现

源码均为本项目自写，Python仅用标准库，C++20使用Apple clang21。
rooted为逐根表，aggregate为四个根矩表。n3/n4两者逐(q,low_b)一致，
rooted还与第40轮独立9!全排列空间账本逐格一致。n5本轮确在MacBook运行，
不是仅有内存计划。保留finished和time回执；该任务exit0且wait/reap，进程已退出。
本机父Agent重新编译并运行两程序n3/n4，输出逐字节匹配保留raw。

以下将PKG设为本目录，OUT设为全新仓外目录，避免覆盖冻结结果：

```sh
mkdir -p "$OUT" /Volumes/data/CodexTemp/cache/clang
env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
 clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pedantic \
 "$PKG/rooted/source/ordered_sector_dp.cpp" -o "$OUT/rooted"
env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
 clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pedantic \
 "$PKG/aggregate/source/ordered_sector_aggregate.cpp" -o "$OUT/aggregate"
"$OUT/rooted" 3 "$OUT/n3-rooted.json"
"$OUT/rooted" 4 "$OUT/n4-rooted.json"
"$OUT/aggregate" 3 "$OUT/n3-aggregate.json"
"$OUT/aggregate" 4 "$OUT/n4-aggregate.json"
python3 "$PKG/rooted/source/audit_ordered_sector.py" "$OUT/n3-rooted.json" \
 --mode n3 --r040-block "$PKG/rooted/references/r040-low-branch-block.stdout" \
 --r040-image "$PKG/rooted/references/r040-image-direct.stdout"
python3 "$PKG/rooted/source/audit_ordered_sector.py" "$OUT/n4-rooted.json" --mode n4
python3 "$PKG/aggregate/source/audit_aggregate.py" "$OUT/n3-aggregate.json" --r069 "$OUT/n3-rooted.json"
python3 "$PKG/aggregate/source/audit_aggregate.py" "$OUT/n4-aggregate.json" --r069 "$OUT/n4-rooted.json"
python3 "$PKG/aggregate/source/audit_aggregate.py" "$PKG/macbook/n5-aggregate-v3.json" \
 --n5-moment "$PKG/references/moment-n5-r001.json" \
 --n5-merge "$PKG/references/first-merge-n5.json"
python3 "$PKG/audit_profiles.py"
```

大内存完整n5重算须在有足够实际可回收物理内存的设备上单worker执行：

```sh
"$OUT/aggregate" 5 "$OUT/n5.json" --plan
"$OUT/aggregate" 5 "$OUT/n5.json" --allow-big-memory
```

n5核心2.625GiB，已记录峰值RSS2,820,407,296字节、墙钟9.92秒、0 swaps。
这里不把小规模重放写成完整n5独立重算。n5逐格的新结果只有这一套全尺寸
聚合生产；旧moment/first-merge提供独立总量交叉，不冒充逐格第二套n5扫描。

原moment-n5-r001.json来源SHA为9f613382e7167fa8dfd94eb9ca71fa938fb7edba1568920c50a4e62db3d3c08e；
旧first-merge输出ee5a100bad13a0ce62676c7d38d9d80c6707ad83b37ef652c9dac4ee241178da，
其历史源码adfa4827b1d3e450193b7351ec2491b21422b7495ebe9076b1a51e43ea742217。
旧源码未改写，来源总量无需旧二进制才能核对。
