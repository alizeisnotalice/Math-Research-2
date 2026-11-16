# MacBook Mathematica GPU：24个40次目标

固定输入seed2026111613、24个首尾1的40次二元多项式，按实根/共轭对分组，交换商第一bit固定0，合计17825792个数值根分配。每批32768，64个半径1.031的DFT采样点；归一为常数项1。全部关键实数矩阵Dot与逐元素运算均检查GPUArrayQ，混合实/复Dot会回退CPU，故实际程序采用实虚部拆分。
Mathematica15.0.0实际113.852秒完成，内核15513显式Exit且PID消失，SSH退出0。最初两次stdin打开方式失败已记录并定向回收；成功方式为WolframKernel -noprompt直接stdin输入代码，不把失败当成功。

stdout压缩账本保留所有批界及2176样本。独立NumPy2.5.4重建全部样本，同时核对根全乘积与目标C绑定：最大根乘积误差7.69e-14、最大样本乘积误差1.12e-7、GPU/CPU系数误差2.34e-6；所有样本违反非负/≤1至少0.34608。只作数值证据，不是全部根分配的严格排除。

复核时先复制本包到CodexTemp并解压stdout.jsonl.gz，然后运行audit.py；其输出写本包目录。准备输入与搜索代码均归档。生产搜索的输入为 data=ImportString[JSON,"RawJSON"] 后执行search.wl，常量inputs.json保证绑定原目标。GPUArray实际设备验证与最终CPU核对均不可省略。
