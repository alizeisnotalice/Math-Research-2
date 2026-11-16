# 短侧内部系数单位传播：q≤20完整条件证书

本包证明在第34轮短侧带符号反设的全部18734个必要模式上，G的每个非零系数都必须为1。不是一般猜想的证明；任意母因子次数由低端/倒数首负界承担。完整规则和最小反例下降见 proofs/035-unit-propagation-descent.md。

输入沿用034-support-patterns/short-G-singletons-q20.json.gz，原始SHA256为06fe63e01e70e621255a7ca33df650cce95d385dc17fa5d8338cabdec5655847。该完整108参数格输入已在第34轮用独立位掩码枚举验证。这里不重复存储。

最终链生成器用五条已证明规则。独立验证器不导入生成器，重建两端状态，逐一核对1,064,638条zero-gap、197,560条zero-unit、137,616条unit、3,401条unit-single-product及263条unit-chain。全部18,734模式覆盖完整，未知0。压缩证书约2MB，展开36.3MB，清单给出两类SHA。早期三规则/四规则扫描分别留下180/14个未知，已被最终链版覆盖，早期中间包留临时目录而不重复归档。

独立Agent以Python3.15.0完成；主线程从归档副本以Python3.13.18重跑28.821秒退出0。负测使用同一验证器，实际拒绝已证零但不整除、可整除但尚未证零、伪造最小支撑、伪造gap见证和删除模式。时间、版本、输入/代码哈希见两个回执。

先将第34轮输入解压到data临时目录，然后执行：

    /Users/zzh/.venvs/main/bin/python verify_unit_chain.py --source /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round034-parent/short-G-singletons-q20.json --chain unit-chain-propagation.json.gz --receipt /Volumes/data/CodexTemp/T2/unfair-01-polynomials/unit-chain-recheck.json

重新发现证书：

    /Users/zzh/.venvs/main/bin/python propagate_units_chain.py --input /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round034-parent/short-G-singletons-q20.json --output /Volumes/data/CodexTemp/T2/unfair-01-polynomials/unit-chain-regenerated.json.gz

输出含运行耗时，重新发现不保证压缩字节相同，应重新逐项验证。只认可验证器的实际运行，不把发现器all_unit标记当证明。
