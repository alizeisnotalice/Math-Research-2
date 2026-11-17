# Round 054：S234 可移植全验证器只读审核

审核对象（父目录，只读）：

- `round054-parent/package/verify_package.py`
- `round054-parent/package/core-manifest.json`
- `round054-parent/package/io-adaptations.json`
- `round054-parent/package/certificate/`

审核时间：2026-11-17（Asia/Shanghai）。本文件只记录审核；没有修改父目录、冻结输入、证书或运行中的任务，也没有启动全量重算。

## 固定版本

| 对象 | SHA-256 |
| --- | --- |
| `verify_package.py` | `b2ad7be4f5874dce2f3a84e986f7c37f4f95203cfe92ca293b93d370a4c0b1cd` |
| `core-manifest.json` | `64d5aeeaa701ad969253ddf797f4779da7e2dbe11a3e362e748ee7478306271c` |
| `io-adaptations.json` | `da95a684a43d6dd84e1d310bf752233da00e488052cd7cd33a7d332e29586721` |
| `certificate/package-manifest.json` | `27ec5c7cca2ec817751da1dd4202d63decd3357d38194d671796c52586674283` |

实际做了以下轻量检查。

- 逐项重算 core manifest 的 212 个文件的字节数和 SHA-256：全部一致，且路径无重复。
- 检查 core manifest：它**包含** verifier 自身、`io-adaptations.json`、所有适配后及 original 冻结代码、58 域输入、树索引和证书；故当前驱动版本不是未绑定的外部脚本。
- AST 解析主驱动与全部 7 个冻结 Python 模块：均成功。
- 检查 `io-adaptations.json` 和逐文件 diff：仅有 `s234_row4.py`、`s234_annulus.py` 的 `RING` 相对路径变换，及下降验证器的 CLI 默认 I/O 路径变换；`s234_row567.py`、`s234_hgq.py`、`s234_ab.py` 字节完全不变。适配后的 `old_ring.py` 也在 manifest 中固定。
- 读取已经完成的 D=61 portable smoke：它从该驱动实际产生 `status=passed`，13 个 survivor 均经 AB/HGQ 双树死亡，共 113 层，最大死亡层 11、最大 live 3。此项是小域运行证据，不替代尚在进行的全量运行。

## 覆盖链审核

主程序的逻辑链闭合如下。

1. 重新枚举固定 58 个 squarefree 域，调用冻结 `R4.one(D)`，并把本次得到的 row4 文件写在新建临时目录；全局只接受固定计数 `ak_pass=61596`、`core_pass=5181392`、`row4_pass=4459150`。
2. 冻结 `R7.one` 只读取这次 row4 临时输出，重新给出行 5--7 和 e/f 的全部状态；全局只接受 `row7_states=5827896`、`state_f_combinations=66621398`。临时 row4/row7 不会作为后续运行的缓存来源。
3. 每个域调用独立下降方向 `RC.one`。它从新生成的 row7 状态重建 base/mother，逐一重算九个正实轴检验与两次 Schur 检验；每个 candidate 的 first-failure byte 都必须等于冻结 gzip ledger，并且重新得到的 survivor 行必须逐序等于冻结 annulus JSON。它不导入历史升幂 root-screen 算法。
4. 对每个 68,593 survivor，主驱动依次重跑 HGQ 树和独立 Fraction AB 树到该行冻结死亡层。它要求两树层级 trace 完全相等（不仅死活相同），再以 `sort_keys=True` 的规范 JSON SHA-256 对照压缩树索引，同时核验 max-live 与最终 state hash。
5. 只在所有 58 个 future 成功完成后才接受全局 4,796,098 bases、54,726,220 mothers、68,593 survivors、621,704 trees layers、`deadmax=71`、`maxlive=7`。单 field 的 `unknown_*`、超时、状态上限、树不一致、byte 不一致或子进程异常均使 `job.result()` 抛出或断言失败，不能生成最终 passed JSON。

因此该脚本不是读取旧 AB/HGQ 成功标志：冻结旧输出只作逐字节/逐 trace 的比较目标，本轮真正重建并执行上述四层必要筛。

## 中断、旧输出与绑定

- `--out` 已存在时立即拒绝；最终结果先写 `.partial` 再原子替换。中断不会留下一个可被误认作正式 `--out` 的 passed 文件。
- 每次运行使用新 `TemporaryDirectory` 生成 row4、row7 与 root audits；不存在可读取的 field cache。checkpoint 仅在某域已经通过本次重算后写入，但该目录在重跑时没有读取分支，故旧 checkpoint 不能造成 false completion。
- 主进程先验证每个 core-manifest 项，结束前再次验证 core manifest 与 driver 自身未在运行期间变化。manifest 对当前的 certificate 输入、冻结 source、适配说明及 verifier 均有明确 byte hash。
- 一个实际发布流程仍应把上表的 `core-manifest.json` SHA 与最终 full receipt 一并纳入 Git/发布清单；这是外层可识别性要求，不是现有驱动的完成漏洞。

## 结论

**接受该 portable full driver 的覆盖和失败语义。** 对固定版本，它足以在一次最终 `status=passed_full_S234_integral_box_prefix_every_root_byte_and_actual_AB_HGQ` receipt 出现时，支持“完整 58 域 integral-box、行 4--7、每个下降 root-screen byte、每个 AB/HGQ prefix layer 均于本次重算通过”的计算性结论。

这仍是该 S234 有限母式/任意有限整数余因子模型的严格程序核验，不应扩大成原始均匀和猜想的一般证明。
