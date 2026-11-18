# 复现第78轮

在本目录用 Python 3.15、选项 -B 运行下列文件。整数 Walsh 系数与二元消元不使用浮点门槛。

| 脚本 | 检验与对应输出 |
|---|---|
| deletion.py，coloops.py | 独立删除与全部四位真值表的仿射层检验，同名 JSON |
| mixed.py | 两个七位谱各自的全部查询核；依赖本目录 audit.py |
| check_weak_folding.py，verify_density.py | 有限域族、R49 与寻址的密度和活动量；weak-folding-results.json、density.json |
| verify_scale.py | 默认最大高度七，输出 scale.json |
| audit.py small，audit.py larger，audit.py mock | 分别对应 small-mini.json、larger-macbook.json、mock-mini.json，主机和耗时会变化 |
| metal_sample.py audit.py ../round077/check_odd_walsh_metal.py | 两台 Apple Silicon 设备各自用 Metal 与直接整数公式对比，metal-mini.json、metal-macbook.json |
| check_lift.py，check_odd_flats.py | 两种多幅共享构造的五代联合卷积；lift.json、odd_flats.json |
| published_seed.py | 公开十位种子的 1024 个完整系数、稠密空间与卷积；依赖 deletion.py |
| gr-verify.py，two_layers.py | 跨层取等、近满核窗口及两层前缀构造；gr-verify.txt、two_layers.json |

MacBook 使用 Python 3.9.6 的项目在实际回执中注明。Apple M6 与 M5 Pro 为不同物理设备。GPU 时间来自 Metal 命令时间戳，整机持续满载率未测量。带系数模板的无界参数结论由轮次证明承担。
