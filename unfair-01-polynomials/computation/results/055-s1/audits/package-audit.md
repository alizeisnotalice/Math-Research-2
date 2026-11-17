# 第 55 轮 S={1} 八次紧凑归档：独立组包回执

组包目录：`package/`。这是索引与证据归档，不是新的全量数学复算，也不代替主 Agent 的 portable-full 检查器。

实际执行：

```sh
python3 /Volumes/data/CodexData/skills/python-math/scripts/run_math.py \
  --workdir /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round055-s1-package/package \
  --log-root /Volumes/data/CodexTemp/T2/python-math \
  /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round055-s1-package/package/verify_archive_index.py
```

运行回执：`/Volumes/data/CodexTemp/T2/python-math/20261117T145835Z-c66d39c3/`；Python 3.15.0，退出码 0。输出严格核对了：1,140 个 core、60,134 条第 4 行记录、67 条第 7 行记录、67 个根筛母式、精确 `H=P\bar P`、根筛失败码、两个原始输出账本、以及 independent/full-v1 的 1,140-core 比对绑定。失败码频数为 `10:4,15:14,33:1,35:2,36:6,37:22,38:12,39:5,40:1`。

最终 `package-manifest.json` SHA-256：

`943c13f2864e6d4283c9741604f7fc4b2b4938d64f76a3abf67e84dbd4260b6c`

检查结果：包内 41 个普通文件（其中 manifest 只是不自列哈希，列出另 40 个文件），0 个符号链接，0 个达到 10 MB 的文件。实际 archive verifier 不含 `/Volumes/data/` 的绝对路径；来源账本中的绝对路径仅用于追溯。

历史披露已写入包内 README：未入包的早期独立 root 草稿误写为 `N*x^2*F`，当时立即不一致，未形成任何通过证书；归档的是修正后的 `F^2+T*x*F+N*x^2` 版本。stage-4 的历史 JSON 整数键 resume 问题同样保留说明，冻结源不修改，后续 portable 验证应 fresh run。
