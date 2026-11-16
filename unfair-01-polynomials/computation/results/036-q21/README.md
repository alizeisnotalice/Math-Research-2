# q=21：独立覆盖与七规则重放

必要模式和非负性结论的前提见第35轮模型。本包不把 unknown 当作可行分解。
完整位掩码枚举独立于发现器 DFS；逐事件验证含七类语义负测。
本轮主线程实际结果见 parent-audit.json：q21有16070模式、16066全部单位；q22有18771模式、18767全部单位；各余4个。后续解析关闭各2个，见第36轮证明。

从本目录复现（临时输出写 data；Python3.13.18，标准库）：

```sh
mkdir -p /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay036/q21
gzip -dc input.json.gz > /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay036/q21/input.json
/Users/zzh/.venvs/main/bin/python verify_q21_two_products.py --source /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay036/q21/input.json --output q21-two-products.json.gz --receipt /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay036/q21/audit.json
```

压缩前后哈希见 input-manifest.json。发现程序保留原路径用于来源追踪；上面的验证命令不依赖原发现目录。q22在主机重启前已落盘完整27.224秒验证回执，原工具session重启后丢失，不伪称重启后收到其退出码。
