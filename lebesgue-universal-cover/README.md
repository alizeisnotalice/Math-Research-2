# Lebesgue 万能覆盖精确面积

研究平面紧凸万能覆盖的最小面积 Λ；目标集合直径至多 1，允许平移、旋转及反射。完成标准为候选覆盖、全局万能性证明及任意覆盖的匹配下界。不得把有限族最优或极限刻画记为原题完成。

材料：STATE.md 记录状态；rounds/ 记录轮次；proofs/ 保存证明；computation/ 保存复现程序与实际结果；references/ 记录来源及本次核验程度。

用户提供的提示词作为研究资料；项目规范及当前授权约束执行。默认三条并行职责为下界、上界、独立审查。远端任务先核验负载，不中止其他课题。

复现：在仓库根目录运行 `python3 -B lebesgue-universal-cover/computation/check_bridge.py`；仅需Python标准库。独立轻量检查：`python3 -B lebesgue-universal-cover/computation/independent_check.py`。

五体放置的严格面积上界：`python3 -B lebesgue-universal-cover/computation/certify_five_body_upper.py 64 4187/5000`。
可选参数依次为每弧段数K（至少64）和严格有理目标；默认值复现第2轮数学证书。
通过已配置SSH别名流式运行完整认证（不在远端落盘）：
`python3 -B lebesgue-universal-cover/computation/run_remote_certificate.py 256 83731/100000 /Volumes/data/CodexTemp/T2/lebesgue-universal-cover/replay-remote.json`。
随后将结果JSON输入 `computation/recheck_certificate.py` 独立复核有理面积、周长及最终比较。
历史证书的源哈希对应各轮提交版本；当前生成器参数化后哈希变化，不意味着历史运行使用了新版。
