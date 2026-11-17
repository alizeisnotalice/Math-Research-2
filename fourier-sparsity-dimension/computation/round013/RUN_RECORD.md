# 第13轮实际计算与复现

日期：2026-11-17。全部归档计算来自本机；MacBook SSH再次在连接阶段超时，本轮没有远端作业，也没有负载读数。未启动MATLAB或GPU上下文。

## 主Agent最终源重跑

capacity_dp.cpp 对给定独立64位真值掩码基，完整枚举码字状态，计算所有支撑子码维数并作锚点覆盖DP。输入含全一字；小定义域函数通过重复坐标嵌入64点，不改变容量。程序检查基独立性、支撑计数为二次幂、全一根状态存在及输出递归有效性。固定使用整数，无浮点容差。

实际编译器：Apple clang version 21.0.0 (clang-2100.4.1.1)，参数 -O3 -std=c++17，显式共享模块缓存 /Volumes/data/CodexTemp/cache/clang。编译产物只在CodexTemp。最终C++源SHA-256：08d26814acc317f1377da6827ca6ec9626f3eb12a2a0e62a9e1ba27e3025ee1d。

主Agent通过 verify_capacity.py 用单个低优先级CPU工作进程串行完成25个实例，总耗时27.373430秒。程序重建基、独立检查全部码字互异、检查分割覆盖和不交、每块指标在码中、并用F₂消元独立核对支撑子码维数和见证得分。每个DP子进程有90秒上限；全部退出0。两个负对照（相关基、缺全一字）均非零退出。

| 实例 | 码维数 | 容量 | 本轮额外完整Walsh核验 |
| --- | ---: | ---: | --- |
| MM k1 | 3 | 4 | 8个相位，原生4点 |
| MM k2 | 6 | 8 | 64个相位，原生16点 |
| MM k3 | 11 | 16 | 2048个相位，原生64点 |
| 常数种子的Q4三角选择器 | 7 | 8 | 128个相位，原生16点 |
| MM k2再作Q2选择器 | 13 | 16 | 8192个相位，原生64点 |
| 20份既有六变量合法码 | 每份15 | 全部16 | 本轮重算容量；bent合法性沿用第11轮完整切片核验 |

共新增10440个小模型的完整整数Walsh相位核验。20份旧码的源输入为 ../round011/period-input.json，SHA-256为445b06403a7479a909f243e43ef059bceab0504ee131f44e5e874cfe4e4f054c；没有把其中旧计时当作本轮计时。最终JSON保存每个实际基、最优分割、局部子码基、时长、退出、源/输入/二进制哈希与软件版本。

分割见证只给容量下界；等于最优值的上界由穷尽DP及证明中的递推承担。25个有限结果本身不证明一般精确倍率定理。

## 第二种完整枚举

enumerate_all_algebras.py 枚举七维Q4三角码的全部29212个线性子空间，筛出115个含单位Boolean子代数。对每个子代数直接枚举乘法条件求最大保留模块，不使用分割DP。最大容量8，共17个最优代数，(代数维数,模块维数)为(1,7)、(2,6)、(3,5)、(4,4)。Python3.15.0实际耗时0.337345秒，正常退出，与DP上界交叉一致。

## 调度、纠错与结束

本机当时有其他研究任务，读取到12逻辑核及较高负载，因此采用单工作进程nice 10串行计算，没有停止其他进程。未宣称本轮本机或MacBook满载。所有最终源运行的子进程均由同步subprocess回收，工具会话返回退出0。

计算Agent早期DP曾错误假设全一字的参数索引，smoke发现后改为实际查找；缺全一或相关基现在显式拒绝。主Agent归档修正源并重新运行全部25例。早期D13和D15计时不计入本报告的27.373430秒；未把临时输出作为最终哈希绑定证据。

理论草稿曾把一般分割单元的共同残余支撑继承列为未知；主Agent用Q±P直接推导纠正，已在正式证明写明。文献备忘中的版本v2与全文已读也经原文复核改成v1和相关页已读。最终材料按实际覆盖记录。

## 复现

从仓库根目录执行，输出放data盘临时目录：

~~~sh
env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -O3 -std=c++17 fourier-sparsity-dimension/computation/round013/capacity_dp.cpp -o /Volumes/data/CodexTemp/T2/a28/r013/capacity_dp_replay
nice -n 10 python3 fourier-sparsity-dimension/computation/round013/verify_capacity.py --binary /Volumes/data/CodexTemp/T2/a28/r013/capacity_dp_replay --period-input fourier-sparsity-dimension/computation/round011/period-input.json --out /Volumes/data/CodexTemp/T2/a28/r013/capacity-replay.json
nice -n 10 python3 fourier-sparsity-dimension/computation/round013/enumerate_all_algebras.py --out /Volumes/data/CodexTemp/T2/a28/r013/algebra-replay.json
~~~

Python程序仅用标准库，主验证要求支持int.bit_count（本次实际3.15.0）。最终文件清单与校验值见manifest.json。
