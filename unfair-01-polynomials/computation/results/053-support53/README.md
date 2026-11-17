# 每侧六点、总次数53的完整认证

本层覆盖2,840,543,082个原始支撑对、1,420,334,032个规范代表；其中1,411,966,297严格排除、8,367,735公平、unknown=0。[证明](../../../proofs/053-six-point-total53.md)保留任意实系数与支撑空隙。结合此前≤52，严格范围推进为每侧≤6点且总次数≤53；九点仍≤34，一般猜想仍开放。

MacBook Python3.12.13、13个CPU工作进程实际独立语义重算全部83,681分片，耗时6652.818秒，最终returncode=0、child_waited=true，父子进程已退出。本机Python3.13.18另以4线程流式读取全部5,420,857,949压缩bytes／78,491,726,768原始bytes，逐片压缩及原始SHA、行数和完整I覆盖全部通过，127.769秒。这是内容绑定审计，不能冒称第二次逐行语义重算。

原始生产是联合53/54任务，中断退出255，保留94,465片。本包只认证实际完整的53切片；旧失败退出与新53成功退出分别原样保存，不把联合任务或54层称为完成。

## 复现与边界

`archive`无损保存83,681份实际语义回执的原始UTF8及SHA，按26个m分包；输入plan、完整原events、payload、旧退出、新53的launch/progress/exit、冻结算法、独立验证器和实际runner均自行封存，不依赖其他轮的路径映射。全部单文件小于10MiB。

```sh
python verify_support53_archive.py --package archive --out /path/to/new-check53.json
```

入口重建全部I计划、核对全部回执及输入绑定、最终会话身份和原始内容审计。它不会再计算14.2亿代表的逐行语义分类；实际逐行分类由已执行的独立验证器承担，SHA a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b，实际runner SHA 7d19ae57220d2d06b0f34c251cbb3decdeae286c3d44ea84bab0da06027a63eb。所有整系数运算及支撑规则均是精确计算，没有浮点容差。

重做全部语义分类时，先从input gzip解压source/plan/payload，按计划在独立临时目录重新生成原始账本，再使用封存的verify_chunk逐片验证。runner中MacBook绝对根目录为运行定位，重新执行时传--root并按其code/plan/run布局放置数据。原5.42GB可再生账本不纳入Git，路径和生成源码已保存；不能只凭本归档验证器的成功重新宣称运行了原始语义计算。

`audits`保留本轮原始字节和语义覆盖回执；`audit_raw53.py`、`audit_semantic53.py`、`build_support53.py`保留实际本机审计/封包过程，其绝对路径对应本轮实跑目录，便携归档验证用上述入口。最终publication-manifest绑定发布文件；哈希不替代数学核验。
