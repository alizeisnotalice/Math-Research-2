# 紧凑支撑枚举与独立复核

两套 C++20 程序覆盖六点总次54、55与九点总次35、36的所有规范支撑对。两者各自产生完整分片，再由 Python 检查连续外层索引、组合数、逐片类别计数与全部未判支撑。两组分片分别为665和1,436片，全部逐记录一致。

归一化后正系数至多为一。只有一种表示的和强制对应的两系数都是一。若一个有多种表示的和包含一对已知单位系数，则该卷积项已用尽预算，额外的正项给出矛盾。无碰撞时各和只有唯一表示，所有正项系数均为一。

生产程序以64位掩码记录相加对角线。并集记录出现过的和，重复位由新平移与先前并集的交集累积。单表示位经反向移位标记两边的单位系数，再检测已知单位对是否落在重复位上。复核程序用位组合后继独立生成支撑，直接在数组中逐对累积表示次数，并用向量词典序判定规范性。两者均枚举共用指数最大公因子为一、较小次数在左、反转规范及同次交换规范的代表。常数因子直接处理。

六点层的生产与复核分别在本机以10进程运行5.67秒和49.79秒。总次54计数与已归档的逐行语义重放完全一致。新增总次55的4,025,316,531个原始支撑对中，规范对2,012,731,962个，单位碰撞排除2,000,704,905个，无碰撞12,027,057个，未判为零。

九点层的生产与复核分别在MacBook以13进程运行7.77秒和166.85秒。总次35有6,197,137,466个原始支撑对、3,098,643,274个规范对、74个未判。总次36有10,641,537,176个原始支撑对、4,895,783,904个规范对、92个未判。全部166个支撑使用163份直接对角支配和3份有理恒等式逐项验过。3份有理证书来自第34轮，按同一支撑重新绑定并重算。各进程均已回收。

从本目录编译后可用标准库 Python 复现。两个输出目录都必须为新路径。按相同方式将总次改为35、36、每侧至多9点及13工人即可重放九点层。

    clang++ -std=c++20 -O3 support_bitset.cpp -o support_bitset
    clang++ -std=c++20 -O3 independent/support_bitset.cpp -o independent/support_bitset
    python3 run_bitset.py --total 54 55 --cap 6 --workers 10 --out /Volumes/data/CodexTemp/T2/a30-six-fresh
    python3 independent/run_bitset.py --total 54 55 --cap 6 --workers 10 --out /Volumes/data/CodexTemp/T2/a30-six-check
    python3 compare_complete.py six54-55-local six54-55-independent --totals 54 55 --cap 6 --out /Volumes/data/CodexTemp/T2/a30-six-compare.json
    python3 compare_complete.py nine35-36-mac-return nine35-36-independent-return --totals 35 36 --cap 9 --out /Volumes/data/CodexTemp/T2/a30-nine-compare.json
    python3 verify_residual.py
