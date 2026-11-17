# 第30轮可携重放

从仓库根目录执行，生成文件全部放到已挂载的 data 磁盘临时目录：

```sh
python3 -B free-boundary-grid-twin-peaks/certificates/r030/replay.py \
  --repo-root . --output /Volumes/data/CodexTemp/T2/a022-r030-replay
```

驱动显式调用同仓库第18、21、28、29轮证书输入，不依赖历史临时目录。它重编
联合DP，顺序重算n3和n5六孔；运行CRT、Pratt、前缀容斥、尾包络、角点插入、
路径核及Fraction检查。原始材料保持只读。编译模块缓存位于
/Volumes/data/CodexTemp/cache/clang-modules。独立路径可用--package-root指定。

本轮实际软件为Python3.15.0（标准库整数/Fraction）和Apple clang21.0.0
(clang-2100.4.1.1)，本机target为arm64-apple-darwin27.2.0。精确任务没有随机种子；
原始Mac模任务的15线程、素数及内存限额保留在raw JSON。已完成的中心四模进程
已退出；后续其他孔队列与两机GPU不在该中心证书范围内。

可携脚本仅更改输入/输出路径接口，原始和归档脚本SHA见IO_PROVENANCE.json。
CRT公开输出仅去除历史路径字符串；输入哈希、Pratt节点、CRT值及整数约束相同。
角点枚举的elapsed_seconds随重放改变，数学字段应保持一致。

首版n3审计曾因JSON字符串与整数直接比较失败；最终审计显式转为整数，未改变DP
或计数逻辑。一般定理由解析推导承担，有限重放不能替代渐近估计。
