# 第55轮：八次单非有理位置1的精确证书

覆盖任意实二次数域、P的唯一非有理位置为1（或反转7）、任意有限带符号整数余因子。其余11个八次支撑轨道和一般猜想仍开放。

核心清单 `core-manifest.json` 固定两套完整前缀程序、两种根筛及预期逐核记录哈希，共7文件；不包含归档材料，避免把索引核验冒充数学重算。

实际全重跑：Python3.13.18，4进程，43.917929666秒，cached=0；1140核心、60134条第四行记录、67条第七行记录，双前缀逐记录及重数相等。67母式全部被精确正轴或外环Schur拒绝，树实际调用0。`full-verification.json`为真实回执，原进程已退出。

复现（标准库即可，输出必须不存在；临时目录放data）：

```sh
python3 verify_package.py --workers 4 --seconds 14400 --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay055-new.json
python3 certificate/verify_archive_index.py
```

`certificate/`无损保存全部必要记录、执行源码及原始来源账本，单文件均小于10 MiB。索引检查只验证归档内容与绑定；第一个命令实际从零重算完整数学覆盖。`audits/output-guards.json`记录拒绝旧输出及`-O`，且原成功回执未改变。

冻结的旧stage4源码存在JSON键/tuple的resume比较不兼容问题，fresh full不读旧cache；历史说明保留在certificate/README.md。一个未发表根筛草稿曾误写范数N项，立即因不一致失败，归档及全部full均使用正确F²+TxF+Nx²。第4行旧摘要文字误称未枚举d，实际源码、完整记录和本次独验均已包含d与Newton4。根筛内环通过/拒绝另有精确fixture，非67实际母式，也非原题候选。
