# 第7项缺口的次数65和66排除证书

[次数方向、支持覆盖与模筛证明](proof-zh.md)。

共覆盖74,448,896个不同二元积，模163的导数 gcd 次数 \(\ge7\) 命中为零。结果允许有理域重因子。既证的低端 \(C_i=1\ (0\le i\le30,\ i\ne7)\) 以及两侧次数至少31是本证书的输入。

本机 Apple clang 21.0.0 与 Python 3.13.18，3 workers。正式全覆盖88.153268秒，含启动和小例的总回执88.940秒，156个工作进程自然退出。峰值单进程1,753,088字节。另实测用含空格异址重编译，从零重建完整支撑文件，复算小盒并独立核完156个原运行分片的哈希。异址验证使用原正式全域账本。

## 重放

在包外取新的输出文件和目录。MacOS 上使用包内入口。

```sh
sh /path/to/package/build.sh
python3 -B /path/to/package/generate_support.py --out /new/path/supports.json
python3 -B /path/to/package/run_sieve.py --input /new/path/supports.json --out /new/path/full-mod163 --workers 3
python3 -B /path/to/package/verify_run.py --input /new/path/supports.json --run /new/path/full-mod163 --compare-stable /path/to/package/frozen-stable.jsonl --out /new/path/verified-coverage.json
```

`frozen-supports.json` 给逐形的允许与必正支持、不相交积立方体及完整计数。`frozen-stable.jsonl` 绑定每片的确定性输出，排除文件路径、PID、内存与耗时。正式按片详细账本留在包外 `full-mod163/shards`，`full-mod163/complete.json` 给汇总。
