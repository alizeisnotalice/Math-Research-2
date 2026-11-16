# 第17轮完成批次：大尺寸权重诊断

所有文件均为完成批次；本目录不包含仍运行的百万粒子工作。n16完整K2用Στ mean_R Zτ归一化，不能除以255；Q为ΣU/ΣZ，公共因子不影响Q。tau2单层不代表完整K2。n16各重复权重高度集中，不能据Q的小数或经验DeltaSE给出认证置信区间。

## 主Agent与独立程序实际复算

运行 `python3 audit_diagnostics.py` 读取本目录5份新原始文件以及已提交r016的两个n12文件，直接用Decimal80解析JSON数值，检查seed、出生秩/重复编号全覆盖、dead零贡献和完整run_end。证书含每层ESS/maxshare及分层delta诊断。所有本批非dead Z均>0；遇到未来活态下溢需要另用log-domain，不可直接套本脚本。

两套源码的本批seed均为 `(0xA022000000000000 + domain + n*1000000 + tau*10000 + replicate) mod 2^64`；domain CLI按十六进制解释。

## 精确输入与复现

n16使用 `../r017-bitset/bitset_smc_checked.swift`，n12使用 `../r016-gpu/production_r013.swift`；编译与输出只写data临时卷：

```sh
WORK=/Volumes/data/CodexTemp/T2/a022-r017-production-repro
mkdir -p "$WORK"
export TMPDIR="$WORK" TMP="$WORK" TEMP="$WORK"
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules
swiftc -O -module-cache-path "$CLANG_MODULE_CACHE_PATH" ../r017-bitset/bitset_smc_checked.swift -o "$WORK/bitset"
swiftc -O -module-cache-path "$CLANG_MODULE_CACHE_PATH" ../r016-gpu/production_r013.swift -o "$WORK/old"
"$WORK/bitset" 16 65536 16 2 256 "$WORK/n16-full.jsonl" 7200000000000000 strat
"$WORK/bitset" 16 65536 256 257 257 "$WORK/n16-one85.jsonl" 8500000000000000 strat
"$WORK/bitset" 16 65536 256 257 257 "$WORK/n16-one86.jsonl" 8600000000000000 strat
"$WORK/bitset" 16 262144 64 2 2 "$WORK/n16-tau2.jsonl" 7300000000000000 strat
"$WORK/old" 12 524288 16 2 144 "$WORK/n12-full.jsonl" full strat 7000000000000000
```

原始JSONL包含每重复GPU/墙钟时间、seed、logZ和完整结束记录。n16完整批1309.88秒；n12加倍粒子批3122.42秒。
`n12-p524288-monitor.csv` 为非独占整机1Hz GPU读数；`resource-exit-audit.json` 保存实际均值99.612%与进程退出检查。
Python3.15、Swift6.4.2、Apple Metal。源码与所有必要输入hash见课题manifest-r017.json。

较小粒子P65536和较大P262144的n16重复都显示明显集中；n12 P524288 Q约.315105，与P262144约.315264接近。这是有限数值对照，不证明稳定极限。尚无混合、权重尾部或有限样本误差认证。
