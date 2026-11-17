# 三次共同因子的独立精确证书

verify_cubics_fraction_bfs.py 只用Python标准库整数/Fraction，实现Sturm枚举、192次有理根二分及逐边有理区间认证。verify_cubic_degree3_portable.wl 独立实现精确根计数、实量词消去与有限状态搜索；未求定的Resolve导致失败。均不读取旧数值搜索或pilot状态。

主线程实际重跑：Python3.13.18，0.1006秒退出0；Wolfram15.0.0，2.8560秒退出0。输出121候选中恰8个必要三次式，四代表状态数180、2、12、2，零态不可达。详细命令：

```sh
/Users/zzh/.venvs/main/bin/python verify_cubics_fraction_bfs.py
/usr/local/bin/wolframscript -file verify_cubic_degree3_portable.wl
```

Python脚本在自身目录写JSON，若重跑应复制到CodexTemp独立目录。Wolfram只打印，允许任意工作目录。父回执包含源hash、stdout、退出码与耗时；independent-final-review记录最初稿审核hash，公开证明追加了全域最小性、一般p与Gauss出口说明。证书只支持明确的条件下降，不解决一般S整数性或原猜想。
