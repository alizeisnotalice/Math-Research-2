# 第24轮实际运行

本机Python3.15.0、Apple M6、macOS27.2、Swift6.4.2；MacBook Apple M5 Pro、macOS27.3、Swift6.4.2。整数/F₂精确运算，零容差。

- expurgate_small_directions.py完整检查15维三角族32768相位，扫描256个商类及各128仿射平移；24坏类被两检查15、224全部排除，留下13维、8192相位、最小重量16。输出含完整基、坏类、保留类。
- 两机相同Swift最终源读取上述最小证书，独立确认15维基、前7个仿射基、两检查秩及核。每机CPU定义式64×64点积与Metal六级FWHT比较524288系数，零分歧，全部相位bent，方向最小重量16。最终输入/源/输出/二进制及命令绑定见receipt。两机均显式指定模块缓存。远端初版漏加固定bent基相位、将方向本身当相位，exit10；最终源已重建z·w并完整重跑。远端前台exit0确认该运行结束，后续SSH进程复查超时，不声称成功完成了事后PID检查。
- explicit_two_checks.py在k2、k3完整枚举原族重量，低重量字分别12、24；k2..5各512个固定种子码内相位全Walsh通过。
- degree_capped_triangle.py完整重量枚举仅四组：(3,1)8192字、(3,2)128字、(4,2)2097152字、(4,3)512字；其余(5,2),(5,3),(6,2)没有全码枚举。全部七组各64个固定种子相位全Walsh检查及显式取等字通过。具体种子、容差、软件及耗时见JSON。

注意重量完整枚举不等于逐字Walsh完整枚举；后者只在明确记录的范围执行。全维合法性、维数、距离由解析证明承担。未打开MATLAB，无持续满载声明。

复验先运行`python3 verify_receipts.py`。Python重跑会更新同目录JSON，建议在副本执行，避免破坏存档输入绑定。Swift复现：

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-module-cache R24ExpurgateMetal.swift -o /Volumes/data/CodexTemp/T2/a28/r024-local/reproduce
/Volumes/data/CodexTemp/T2/a28/r024-local/reproduce expurgate_small_directions.json /Volumes/data/CodexTemp/T2/a28/r024-local/reproduce.json
```
