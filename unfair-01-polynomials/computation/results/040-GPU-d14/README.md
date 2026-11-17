# 第40轮：d14恢复完成与保存行复核

实际参数 m=60,n=100,d=14,starts=32768,steps=12000,restart=2000,margin=.001,seed=2026111615。由第39轮记录的step7000/epoch3真实state恢复，最后完成12000步；wrapper退出0并wait回收，主线程ps确认25270/25272均不存在。release标记只作内核内部记录，实际退出另由queue及ps核验。

`run/`保存实际config、五份restart、最终checkpoint/result、CPU前缀审计、恢复来源和release。全部96行JSON系数独立以Fraction重算；主线程从本归档路径再次完整重建前缀池并重算96行，得到与独立审计逐字段相同的`parent-fraction.json`。保存的小数按其精确十进制有理数理解，最小残差约0.001765378750860691，均不是二元乘积的精确反例；不能据此排除未保存或未搜索对象。

真实恢复含step7000旧记录及同一步的`resume_initial_evaluation`，history不是严格递增。绑定审计明确接受且单列这唯一有来源重复，其余下降或重复均拒绝。没有将它误记为第二次优化推进；未来重启随机流也不声称与断电前逐比特一致。

冻结worker和prefix源码在`../039-d14-recovery/`。复现不需GPU：

```sh
python verify_prefix_fraction.py --run-dir run \
  --worker ../039-d14-recovery/resume_gpu_factor_search_prefix.py \
  --prefix-module ../039-d14-recovery/frozen-worker/fractional_prefixes.py \
  --record parent-fraction.json
```

`audit_d14_resumed_binding.py`另核原始恢复状态、真实Adam终态及来源哈希；其需要本地原始state文件，不伪称公开JSON等同51MB的tensor状态。该状态与原始恢复链保留于项目.local归档，目标盘备份尚未核验。只上传可复验的最小数值记录，不上传大tensor缓存。
