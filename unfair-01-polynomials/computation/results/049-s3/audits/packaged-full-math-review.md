# 第 49 轮封装 full 验证补充审计

本补充只审计第 49 轮 S3 证书封装后的实际 `full` 回执及新增发布材料，不重新执行其
全量计算。

## 已核验通过

* 回执
  `/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round049-parent/packaged-full-verification.json`
  的 SHA-256 为
  `04e8872757383294a08b9797185099c9ec0c5369ec40bfaca487952acd8e88b6`。
* 它记录 `mode=full`、`screen_replayed=true`、`workers=4`、Python 3.13.18、实际
  `verify_package.py --mode full --workers 4` 命令，以及
  `elapsed_seconds=355.0039487499962`。
* 其严格汇总为：669 个 TN、27,909,836 个 Newton 母式、90,092 个幸存母式树、854,628
  个树层、最大死亡层 53、最大活态 27，状态为 passed。
* 回执中的验证器 SHA
  `ae85219ef7bfda8b06fddd7975b7e0fe8f3976db4ca1fdef01407b85a9ac6f67`
  与读到的并行版 `verify_package.py` 一致；manifest SHA 也记录为
  `31403d79ebea8495bfb9c8d8d67ac03a31e69329a3656ed28fae354f38025983`。
* 已用 `ps -p 91639` 核验启动 full 的 PID 不再存在。并行器按固定 669 索引提交 futures，
  每项经 `future.result()` 后才可能进入最终断言及写出 passed 回执，故 worker 异常不能被
  混记为通过。
* 新 README、历史完整驱动和历史独立运行回执实际存在；其 SHA-256 分别是
  `fe1579e771ac4ae614aecb0560fcd42686a9a9ceac8d64da87b78a6bc37141de`、
  `c1f6fa27b30f77a2be943868100a553d6b2955e285051d575db5c8ea9e860fd9`、
  `1ac65f79802cc26adab2b7757bd6ab6c40816c4335711d7d338c2cbb9d4255e7`。

## 尚需一项发布索引修正

当前 README 声称存在 `publication-manifest.json`，用以绑定 README、验证入口、运行回执和
数学审查；但本次 `rg --files` 与直接读取均确认包根目录没有该文件，且原 `manifest.json`
未列这些新增非核心发布材料。故核心数学证书和 full 回执均已通过，但发布索引的该一句
尚不自洽。

发布前应二选一：建立并列出所称 `publication-manifest.json`，逐项哈希绑定这些材料；或删改
README 的该项声明，改为实际存在的可核验索引。完成后，本轮先前提出的 README/独立全筛
溯源条件即可关闭。
