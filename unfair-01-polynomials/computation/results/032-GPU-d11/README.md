# d11完整批次与保存样本独立审计

实际参数m48,n88,d11、starts32768、steps12000、restart2000、margin.05、seed2026111609。六个完整cohort共196608个轨迹初值，每个cohort覆盖3181个必要前缀，各前缀分配10或11行；独立前缀枚举为59049原始编码、5131低端合法、3181筛后编码。

实际运行4940.674秒。session74904退出0、父44394/子44395均已消失并回收；`d11-exit.json`保存实际检查。release中的103023360字节是退出前的tensor记录，不能解释为退出后仍有常驻内存。随后才启动独立d12批次，运行快照另存。

主线程从本归档副本实际复跑下列两个标准库程序，退出0：

```sh
python3 audit_prefix_coverage.py --run-dir run --worker worker/gpu_factor_search_prefix.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/d11-coverage-recheck.json
python3 verify_prefix_fraction.py --run-dir run --worker worker/gpu_factor_search_prefix.py --prefix-module worker/fractional_prefixes.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/d11-fraction-recheck.json
```

它们重新生成前缀池，并对五次重启各16行及最终16行共96个唯一保存样本以Fraction重新卷积。最低保存残差为0.05019746720790863，所以这些保存小数中没有精确二元分解。浮点搜索和未保存轨迹不因此被排除；margin.05也未覆盖接近端点的系数。

审计代码复用第30轮已核验的通用算法；覆盖程序仅更新描述文字，Fraction验证器逐字节相同。运行模型PyTorch2.15.0、NumPy2.5.4、MPS float32；本次父标准库审计Python3.15.0。`run/`只归档必要JSON，可再生的state.pt保留在临时目录。checkpoint未保存全RNG，不承诺逐位续跑。
