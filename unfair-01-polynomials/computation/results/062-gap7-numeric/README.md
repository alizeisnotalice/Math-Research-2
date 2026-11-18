# 本机前缀预算 MPS 批次

次数 \((24,40)\)，种子 2026111811，16,384 个并行初值，32,000 步。源码沿用[第61轮预算参数化](../061-gap7-numeric/search_gap7_mps.py)。全部 256 个保存样本用整数分母卷积验证，记录的最小最大残差约 0.24504。已保存结果、释放缓存并退出。

    python3 -B ../060-gap7-numeric/verify_saved_samples.py --out /Volumes/data/CodexTemp/T2/a30-gap7-m24n40.json m24n40

另一批次数 \((28,44)\) 把两个首项饱和所迫的零位写入 GPU 参数化，同时优化平方均值和最大平方残差。种子 2026111812，16,384 个并行初值和 36,000 步。全部 288 个保存样本有理复验，最小记录最大残差约 0.21557。已保存并退出。

    python3 -B ../060-gap7-numeric/verify_saved_samples.py --out /Volumes/data/CodexTemp/T2/a30-gap7-m28n44.json m28n44
