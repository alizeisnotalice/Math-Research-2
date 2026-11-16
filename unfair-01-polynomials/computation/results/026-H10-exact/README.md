# H10(2,3) 完整精确证书

本包仅闭合 H10(x²)H10(x³) 固定目标。8310 个必要根布局中8309个有严格负系数，1个为标准分解。覆盖证明在 ../../../proofs/026-H10-root-layout-coverage.md。

从本目录运行：
```sh
python3 verify_H10_exact.py --input H10-a2b3-inputs.json --balls H10-a2b3-root-balls.json --certificate H10-a2b3-certificate.json --galois H10-mother-certificate.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H10-recheck.json
python3 test_H10_checker_guards.py --input H10-a2b3-inputs.json --balls H10-a2b3-root-balls.json --certificate H10-a2b3-certificate.json --galois H10-mother-certificate.json --scratch /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H10-negative-recheck
```

验证器只用整数与Fraction，不依赖options缓存或搜索器。两份完整验证回执及主线程14项语义篡改测试均为实际执行结果。原独立11项回执保留，当前测试脚本为主线程修订的14项版本；重复母根测试同步Rouché数据以实际抵达分离检查。

生成脚本位于 generation/，其root_ball_common.py与验证器的精简独立helper刻意分开；原始103MiB分数区间可再生中间件不入库。证书 compact v1 的边界是有理二进区间，checker重新计算每个指定系数，不信任存储符号。环境：Python3.13.18；生成使用mpmath，验证无第三方依赖。SHA256SUMS覆盖本包最终文件（不含自身）。
