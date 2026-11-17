# 第49轮本机深窗口96样本

本机 MPS 实际完成 m80/n120、d15、B8192、36,000步、每6000步换窗、seed2026111705。六个窗口覆盖既有46,011个前缀，各保存16行，共96行。原进程76790已退出，release记录完整结束；没有把预估路径当作GPU实跑。

本包含冻结worker、base、前缀工具和计划、实际config/结果/窗口样本及父独立Fraction验证器。保存的十进制数被视作精确有理数重算全部卷积，96行均不是精确反例，最小最大残差约0.0011671929950154425；具体有理残差见回执。这只核验保存行，不核验整个训练轨迹，也不证明任何有限范围排除。前缀池的完整性沿用第47轮证书，不重新扩张其模型。

在本目录，用Python标准库重放，无需GPU、PyTorch或原临时绝对路径：

```sh
python parent-generic-audit.py --run . --verifier fraction-verifier-from-r040.py --worker worker.py --m 80 --n 120 --seed 2026111705 --margin 0.001 --prefix-plan prefix-plan.json --base-source frozen-base-worker.py --out /path/to/deep96-audit.json
```

config保留原运行路径用于溯源，命令显式覆盖依赖路径。parent-m80-deep96-audit.json为原本机父核验，portable-deep96-audit.json为独立目录封装重放，publication-audit.json为最终公开目录的父再次重放。源/输入SHA、Python实际版本和逐行残差均保存于回执。GPU worker需实际PyTorch MPS才能再次训练，不能将上述精确验证命令描述为GPU训练。
