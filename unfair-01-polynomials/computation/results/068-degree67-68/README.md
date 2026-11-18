# 第7项缺口的总次数67和68排除证书

- [次数方向、完整有限覆盖和算术筛](proof-zh.md)
- [任意次数的支撑传播定理](support-theorem-zh.md)

使用已证的 \(C_0,\ldots,C_{30}\)（\(C_7=0\)，其余为1）、\(I_{31}>0\) 及两侧最低次数32，生成81个支持分支并精确取并。总次数67覆盖23,592,960个不同积，68覆盖125,304,832个不同积。模163导数 gcd 次数达到7的命中为零，结论允许有理重因子。

Apple clang 21.0.0 与 Python 3.13.18，3个本机工作进程、nice+4。全域阶段耗时180.830352秒，含小例的外层回执182.208秒，单进程峰值1,753,088字节。293个工作进程自然退出。另在含空格异址重编译，重建全部支撑、通用传播回归和小盒，并复核正式全域账本。

## 重放

在包外选尚未存在的路径。编译共用数据盘模块缓存，Python 需要启用断言。

```sh
sh /path/to/package/build.sh
python3 -B /path/to/package/generate_support.py --out /new/path/supports.json
python3 -B /path/to/package/verify_general.py --input /new/path/supports.json --out /new/path/general.json
nice -n 4 python3 -B /path/to/package/run_sieve.py --input /new/path/supports.json --out /new/path/full-mod163 --workers 3
python3 -B /path/to/package/verify_run.py --input /new/path/supports.json --run /new/path/full-mod163 --compare-stable /path/to/package/frozen-stable.jsonl --out /new/path/verified-coverage.json
```

`frozen-supports.json` 保存各方向逐形支持、固定位及不相交立方体；`frozen-general.json` 保存通用规则的81形复算及小抽象输入；`frozen-stable.jsonl` 保存每片的确定性账本。正式按片完整输出在包外 `full-mod163/shards`，运行汇总为 `full-mod163/complete.json`。

MacBook 完整复跑 293 片，模 163 零命中，各片稳定散列逐项相同。回执为 mac-complete.json 和 mac-verification.json。主进程另核对远端生成的支持表与冻结文件字节一致。
