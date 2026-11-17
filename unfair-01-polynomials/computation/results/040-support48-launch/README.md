# 第40轮 MacBook 47／48 生产批次

这是已实际启动、尚未独立数学认证的生产批次。参数为两侧支撑至多6点，总次数47或48；47个次数对、100005个I分片，15个fork工作进程。每批提交上限14400秒，到时停止提交、排空并join；成功输出才记完整生产，超时或失败仍为部分未知。

`plan/support48-plan.json.gz` 是实际完整输入清单的确定性 gzip 副本；原文 SHA-256 为 `ba259776db7eb704c7411c9d5c66739736aede2d9e091ab62ed52b448e61a04d`。`support48_prepare.py` 从独立组合生成清单。`test_support48_plan.py` 另行重建输入覆盖、运行两项真实chunk并拒绝重复/错误I/越界输入。

本目录的 `remote_support48.py`、`support48_bridge.py` 是实际生产源码。桥中机器别名、已验证地址、解释器和输出路径属于本机配置；移植时须使用自己的已核验SSH配置。归档不含凭据。原本题GPU任务已结束的资源另行回收，其他课题的GPU任务不受本CPU批次影响。

在临时副本解压后可复核：

```sh
gzip -dk plan/support48-plan.json.gz
python test_support48_plan.py
```

实际启动与收到的中间计数见 `launch-and-progress.json`。其中46记录是42553个缺片已生产完成、旧1998片加新片共44551；远端旧父子PID已核验消失、池已join。46独立逐项重验仍在运行，本轮不登记≤46定理。47/48的中间计数同样不是完整覆盖，也不代表严格排除。
