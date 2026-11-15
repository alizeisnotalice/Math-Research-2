# n=7 beta=0 分轨可复现实跑认证

本次运行的目的，是建立一个有运行前源码—二进制哈希对应的 n=7 精确复现记录，并逐批退出以释放内存；它不改变 beta 递推的数学结论，也不构成渐近结论。

## 冻结输入与执行方式

运行前冻结源码 `beta_zero_layered_orbit_certified.cpp`：

```text
e0c769ad793d63391dc3e853b7170680afa04dfde4b53c1c656f9d1fd875c3b5
```

以 Apple clang 21.0.0、C++17、`-O3 -pthread` 编译为唯一命名二进制
`beta_zero_layered_orbit_certified`，运行前 SHA-256：

```text
945b1cb671f8c04bc36c5e267de79f1e3c73783644f4b2e3ef5a7ae88d622cb3
```

串行驱动 `run_n7_beta_orbits_certified.py`（SHA-256
`b2538d78aec76b2c14944cb299dec544c46b7ce0ca8fb0ae4e0b3021ecb285ee`）
依次运行 `--orbit 0` 至 `--orbit 11`。每一进程内部使用 12 workers；驱动等待该进程完整退出、检查其 JSON 后才启动下一轨。运行前/后源和二进制哈希均相同，完整命令、时间戳、退出码、JSON、stderr 和 `/usr/bin/time -l` 输出见 `beta-zero-layered-orbit-certified-n7/manifest.json` 及同目录的每轨文件。

## 结果与独立核验

12 个轨道均 exit code 0，均通过退出后 JSON 自检。独立标准库验证器
`verify_beta_zero_layered.py`（运行前 SHA-256
`87c916e2f7e36447b2204b7d352f1fc4e1f13ce9a9c4e2cabee79c9bff6e4d8e`）
检查：

- 共 12 个无向 D4 轨道及有向边覆盖 168；
- 每轨资源文件包含 `maximum resident set size` 和 `peak memory footprint`；
- 逐轨加权和严格等于
  \[
  B_0(7)=131529421404957780122938007473728844832061946591752192;
  \]
- 它与历史完整运行 JSON 的总数完全相同，且 n=3/5/6 的既有独立交叉值也都通过。

实际独立核验证书为
`beta-zero-layered-orbit-certified-verification.json`，SHA-256：

```text
baff028dbe0169e3cd8d22b764c53a6d34212178f9b8261b197ead23d396795a
```

该结论仍只确定有限 \(n=7\) 的 \(\beta=0\) 与 \(\beta>0\) 支持质量，不能给出 \(\mathbb E_2\beta\)、\(R_n\) 或题设极限的渐近断言。

## 每轨实际资源

| orbit | max RSS (bytes) | peak footprint (bytes) |
| --: | --: | --: |
| 0 | 4,445,569,024 | 8,364,300,576 |
| 1 | 6,517,899,264 | 8,360,319,360 |
| 2 | 6,048,743,424 | 6,051,173,128 |
| 3 | 6,102,319,104 | 8,779,192,968 |
| 4 | 5,911,789,568 | 6,944,265,568 |
| 5 | 6,491,815,936 | 6,805,394,712 |
| 6 | 5,645,287,424 | 5,735,731,712 |
| 7 | 5,348,950,016 | 5,405,937,960 |
| 8 | 5,148,065,792 | 6,127,981,392 |
| 9 | 5,311,578,112 | 6,167,024,440 |
| 10 | 5,148,409,856 | 6,280,860,616 |
| 11 | 4,979,556,352 | 6,203,036,496 |

内存峰值在单轨之间不累积；这是按独立进程顺序执行而非此前单进程跨轨路径的实际记录。认证结束后已检查，没有 `beta_zero_layered_orbit_certified`、驱动或验证 Python 进程残留。

## 主要文件

- 输入与运行前哈希：`beta-zero-layered-orbit-certified-input.json`。
- 冻结源码、二进制：`beta_zero_layered_orbit_certified.cpp`、`beta_zero_layered_orbit_certified`。
- 顺序驱动：`run_n7_beta_orbits_certified.py`。
- 逐轨输入/输出/资源/退出记录：`beta-zero-layered-orbit-certified-n7/`。
- 聚合验证器与证书：`verify_beta_zero_layered.py`、`beta-zero-layered-orbit-certified-verification.json`。

历史单进程运行的二进制哈希缺口仍保留在较早报告中；本认证不倒填或声称恢复该历史哈希。
