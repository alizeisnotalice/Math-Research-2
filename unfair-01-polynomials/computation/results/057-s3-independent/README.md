# 位置3分支的独立根筛核验

正轴部分以齐次整数评价 \(2q^8P(p/q)\)。其最小嵌入等于整数线性式减 \(p^3q^5\sqrt\Delta\)。故单点给出一条关于 \(d,e\) 的严格条件，等价于该整数线性式至少为 \(\lfloor p^3q^5\sqrt\Delta\rfloor+1\)。对原压缩组按每个 \(d\) 交 \(e\) 的整数区间，全量重算 13,209 核心、1,915,327 组及 1,523,605,945 个母式，得正轴保留864,141,765项。主程序另以 Horner 整数式和平方比较复核 5,216,058 项，全部相同。三个试跑核心的3,788,628项全部逐项交叉核对。

根环另由 MacBook 13 进程应用二进制定点整数区间。二次根区间通过整数平方二分取下界和上界。每级 Schur 递推用有符号128位四端点乘积和区间减法，再作向外舍入的二次幂缩放。首项下界严格大于常数项绝对值上界才继续。首项上界不大于常数项绝对值下界时可严格拒绝。其余情况计未判并保留。本批未判0，四种首失败计数均与主生产一致，869,396个保留母式完全相同。MacBook执行3.653秒，主管等待13个子进程退出。

独立整数范数校验从每个根筛通过项重建 \(H=P\bar P\)，分别以升降幂两向严格检查根环。869,396项均通过，两侧根环、原Newton尾区间成员资格和整数范数指纹均已保存，实际70.087秒。

复现定点区间计算可使用主证书生成的规范文本尾组、split_s3_groups.py 及以下命令。从新目录启动，先复制 run.py 与源码，编译后再运行，保留原只读回执。本归档的13路通过表按原字节压缩，复核器会重新解压、核SHA、合并并逐项比较。

    python3 split_s3_groups.py /Volumes/data/CodexTemp/T2/path-to-fresh/full-groups.txt /Volumes/data/CodexTemp/T2/a30-s3-integer-replay
    clang++ -std=c++20 -O3 s3_integer_enclosure.cpp -o /Volumes/data/CodexTemp/T2/a30-s3-integer-replay/s3_integer_enclosure
    python3 verify_s3_screen_comparison.py --production /Volumes/data/CodexTemp/T2/path-to-fresh/full-pass.jsonl --return-dir integer-remote --manifest integer-remote/manifest.json --source s3_integer_enclosure.cpp
