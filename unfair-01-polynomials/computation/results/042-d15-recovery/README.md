# d15 真实优化器状态恢复

本适配器已经实际启动，PID93524；旧PID74343为零RSS、无文件描述符的退出残留，确认不再执行后才启动。state.pt的可靠边界为1500，旧checkpoint中1600/1700被丢弃；不是从较晚浮点摘要伪造优化器状态。源worker/plan/order/Adam检查见实际provenance和CPU回执。后续已真实完成并保存resumed_progress；本包progress-snapshot仅为运行中时点，不代表批次完成。

SOURCE常量绑定第41轮归档的原worker字节，运行时其本机路径需存在。生产复现参数：m64,n104,starts49152,steps12000,fractional-index15,margin.001,restart-every2000,seed2026111701,seconds14400；原plan与输入在041-GPU-d15-input。给定 --resume-run、--resume-state、--prefix-plan 和新 --out，先用 --validate-resume-only 验证，再启动MPS。具体参数名可由 --help 检查；不向原run写入。

实际短程CPU负测拒绝错误张量形状和Adam学习率。fsync与原子rename保存状态，期限到达记partial_time_limit，异常记failed_exception；未保存完整RNG，未来随机重启不承诺逐位连续。本包不含80MB state，必要状态另存本地.local并核对hash；这不等于目标盘备份。
