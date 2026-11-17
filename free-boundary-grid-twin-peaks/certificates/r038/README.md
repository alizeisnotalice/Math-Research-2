# 第38轮证书与复现

prefix-occupation包含精确占用恒等式、逐行反例、两种几何权独立核验器及有限结果。finite-words为四个完成搜索链，不代表整批完成。eleven-prefix包含完整613分片账本、独立覆盖检查与全N解析提升；冻结计算源沿第37轮。

```sh
python3 -B free-boundary-grid-twin-peaks/certificates/r038/replay.py \
 --output /Volumes/data/CodexTemp/T2/a022-r038-your-replay
```

输出目录必须全新且位于data卷临时目录。需要Python标准库和clang++；共享编译缓存显式位于data卷。重放10项有限检查，不重复完整N11枚举，也不启动生产搜索。
