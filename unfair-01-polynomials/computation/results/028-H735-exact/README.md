# H7(3,5) 固定目标证书

证明见 `../../../proofs/028-H735-exact-theorem.md`。全部119062必要布局中119061严格负，唯一标准mask50792199。一般猜想未解决。

主线程在归档副本实际完成下列全量认证（Python3.13.18，标准库整数/Fraction；7.47秒，退出0）：

```sh
python3 verify_h735_ball_certificate.py --balls H7-a3b5-root-balls.json \
  --masks canonical-masks.json --certificate certificate.json \
  --galois H7-mother-certificate.json \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H735-recheck.json
```

重新从所有29位群选择验证根盘、必要布局及外部独立账本：

```sh
python3 verify_frontend.py --balls H7-a3b5-root-balls.json \
  --source-ledger H7-a3b5-canonical-masks.jsonl.gz \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H735-front-recheck.json
```

`test_ball_kernel.py`、`test_compact_schema.py`、`test_certificate_guards.py`、`test_frontend_semantic_tampers.py` 已从本目录实际运行；测试输入路径改为本包本地文件，完整验证器、前端、核函数保持原始字节。前端语义测试使用 data 内临时目录并自动删除。

`certificate.json` 仅把 side/index 当作证书建议；验证器独立重建96位整数球，不导入搜索器。源精度128位与验证精度96位分开记录。12项额外 Fraction 抽检不是全量认证的替代。最小严格负余量见两个完整验证回执。

若重跑搜索，必须指定新空临时目录：

```sh
python3 search_compact.py --balls H7-a3b5-root-balls.json \
  --parent-masks canonical-masks.json \
  --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H735-new-search
```

根中心生成器使用mpmath110位，再由有理Rouché验证。重生成应复制脚本和helper到独立临时目录运行，避免覆盖证书。`build_h7_a3b5_layout.py` 保存原始独立枚举实现及来源绑定，原始版本仍含本机母证书定位；严格复现认证不依赖再次运行该生成器。13MiB未压缩原始账本未提交，gzip解压哈希为928a80ae9b0b405db118bd5d01945f7c9c6402cf5c644488ff2edccbdc6fcb0c。

`canonical-masks.json` 文件含末尾换行，SHA为0f89d66edc9656f0610b6e647430fbc23d073e6ed847f439d248ea6d4d96278d；证书绑定的是无末尾换行的规范编码，SHA为9ddefd4f202e2fe9bf36be38d685096bc5abb9a298779f7b74a44703d62c87d2。二者不混用。`source-map.json` 保留原文件SHA与测试路径适配说明。
