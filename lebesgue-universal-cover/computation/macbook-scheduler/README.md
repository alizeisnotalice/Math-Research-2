# Lebesgue 的 MacBook 常驻队列

本目录含一个 Python 调度器和一个重复使用 Metal 设备、命令队列及共享缓冲区的七体搜索后端。调度器持有项目根目录中的排他锁。通过 detach 启动时，进程在独立会话中继续工作。自定义 Metal 搜索共用同一个串行 GPU 名额。七体后端空闲时保持驻留，自定义 GPU 作业运行期间它不提交命令。

## 状态与恢复

在远程项目独立目录将以下路径统一置于 ROOT 下。

| 路径 | 内容 |
| --- | --- |
| scheduler.pid、scheduler.lock | 调度进程号和单实例锁 |
| ledger.json | 稳定任务键、领取状态、活动 PID、结果散列、遗留作业完成状态 |
| heartbeat.json | 任务队列、当前压力、CPU预算、GPU进程、是否需补队列 |
| inbox/*.json | 只追加的批次文件 |
| runs/任务键/attempt-NNN | 命令输入、stdout.log、stderr.log、worker.json、status.json、验证过的结果 |
| gpu-spool/ | 持久七体队列、设备日志、压力提示 |
| STOP | 停止领取新任务的标记 |

账本、心跳和完成状态采用同文件系统的临时文件及原子改名。只由调度器写账本。子进程在独立会话内完成尝试并留下状态，调度进程重启后先核对已有尝试。进程退出且没有可靠回执的任务保留为 interrupted。失败或中断可以显式重领，已完成键保持唯一。

SSH 恢复后先读取 heartbeat.json 和 ledger.json，再核 PID 和文件。保持同一个 SSH 控制套接字，已经有一条远程连接且仍能使用时，后续通过它发出状态读取和上传。日常换轮通过 inbox 追加新批次。达到阶段目标或用户暂停时创建 ROOT/STOP。调度器停止领取，允许当前尝试保存回执后退出。支持协作保存的自定义任务可设 checkpoint=true 并监测环境变量 LEBESGUE_STOP_FILE。它指向该运行目录中的请求文件。旧外部作业由原拥有者控制。

恢复前核对旧进程已经退出并移除 STOP，然后执行

    "$PYTHON" -B "$ROOT/scheduler.py" detach "$ROOT"

带有同一锁的第二个 serve 实例以73退出。MacBook 关机重开后可按上式从原账本恢复。无需修改系统登入设置。

## 提交格式

先上传脚本和不可变输入，再把清单写在 inbox 内同目录临时文件。写完后原子改名为一个新的 JSON 文件名。每项的稳定键为 stage、version、input_sha256、scope 四元组 JSON 数组的 SHA-256。改变算法、输入或范围应更新相应字段。可用下列指令从单项 JSON 预先计算依赖键。

    "$PYTHON" -B "$ROOT/scheduler.py" key "$TASK_JSON"

提交一个自定义 Metal 作业的格式

    {
      "tasks": [{
        "kind": "gpu",
        "stage": "counterexample-scout",
        "version": "algorithm-v1",
        "input_sha256": "<真实输入的64位小写 SHA-256>",
        "scope": "明确宿主范围和手性",
        "priority": 50,
        "cpu": 1,
        "cwd": "/absolute/project/workspace",
        "argv": ["/bin/zsh", "/absolute/project/workspace/run.sh", "{run}/result.json"],
        "result": {"path": "result.json", "type": "json", "equals": {"complete": true}}
      }]
    }

CPU 任务使用 kind=cpu，并在 cpu 中填入验证器真实使用的线程数。调度按核预算分配。数学库的嵌套线程默认一核。argv 按字面参数执行，支持 {run}、{root}、{threads} 替换。由上传者对自定义命令及输入负责。若任务从标准输出生成完整回执，可设 stdout=true，调度器在退出后原子发布。输入若由七体后端直接读取，调度器会核对文件 SHA-256。

路径 result.path 相对本次 attempt。type=json 按 equals 检查指定字段，也支持以点号访问嵌套字段。type=jsonl 可用 rows 指定完整行数，用 row_keys 检查逐行字段。退出码和结果内容都通过才标记 complete。每个 after 键都处于 complete 后，依赖任务才有资格领取。回执完整性由调度检查，面积与覆盖等数学结论仍由各任务的专用验证器判定。失败任务的重领文件格式为

    {"tasks": [], "retry": ["<已有失败或中断的任务键>"]}

七体原生任务把 argv 换为 driver=seven、input（绝对文件路径）、epochs、seed，result.path 固定为 result.json。每份输入是原算法规定的14个姿态坐标及82个点坐标，用空白隔开。GPU 回执包含 Metal 设备、实际命令缓冲数、已计分姿态数和设备用时。它给出浮点候选排队信息。seven_followups.py 用 after 绑定其结果，自动把已完成的七体候选交给 Arb 固定摆放核验。

## 资源和连续续接

在 config.json 配置 cpu_budget 和 metal_source。每次心跳核对作业及实际进程。Linux/macOS 的常规数学子进程由进程关系识别。MacBook 的内存压力由 memory_pressure、交换量变化和系统热限制合并。系统升压时减半新派发额度，严重时暂停新派发并向支持检查点的高占用研究任务写请求。既有外部作业只观察。原生 Metal 每个提交块的子步数量随压力调小。CPU 严格任务优先于低优先级候选整理，有 GPU 工作待领时预留 GPU 供料核。

external.json 的 jobs 项用 id、pid、command_contains、cpu、gpu 描述已有作业。完成契约可指定绝对路径、内容条件及 exit_path。调度在旧父进程或已观察的旧子进程活动时保留其额度。旧 GPU 作业占用未清空前，统一队列不领取新 GPU 作业。缺少完成契约的旧退出记录只标为 ended_unverified。

七体备用搜索冻结已完成的真实原叶结果，按分片轮流选不同构形。每波八个独占 GPU 子批，前六个提交结果通过后，调度执行一个短 CPU 任务生成下一波清单。每个 GPU 子批另派生一次 Arb 核验。任务范围耗尽时设置的心跳补队列提示通知研究工位补新输入。上界和反例可用更高 priority 插入同一队列。已开始的有效批次先完成。

## 初次布置与验证

用已经核实的远程 Python 环境，保持该项目原 seven.metal 不变。复用已配置的 Clang 模块缓存。

    xcrun clang++ -O2 -std=c++17 -fobjc-arc -framework Foundation -framework Metal \
      -fmodules-cache-path="$CACHE" "$ROOT/resident_gpu.mm" -o "$ROOT/resident-gpu"

bootstrap_scheduler.py 绑定本次交接的旧七体和0.83741原树任务。它适用于具备这些历史任务资料的同一研究目录。随后可用七体已有定点验证脚本运行 seven_followups.py 的 produce，和七体队列入口 seven_wave.py ROOT 0。其他研究线路各自上传清单即可。

在隔离的临时目录可测试排他锁、错误回执、依赖门及 STOP。

    "$PYTHON" -B smoke_scheduler.py "$SCRATCH" 64

此测试使用小型 CPU 文件作业。实际 MacBook 交接另验了两路旧 GPU 任务、一个旧候选 CPU 搜索和三个0.83741验证作业的完整回执。随后新 Metal 与 Arb 任务均经同一个账本，由已分离会话的服务继续派发。
