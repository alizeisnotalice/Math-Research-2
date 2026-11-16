# 九点33／34层：完整覆盖及全部精确余例

本包支持[九点总次至多34定理](../../../proofs/033-nine-point-total34.md)。17片、33个次数对完整；6,622,858,623原始支撑对约化为3,138,096,530规范代表，46个余例现全部有精确证书。

九点账本原始压缩文件保留于 /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round020-fullcpu/output/nine34-m01 至 nine34-m17。check_nine34_receipts.py 实际完整解压195,950,258,308字节，逐片哈希与既有独立逐项验证回执一致。轻量补充验证器另核对源程序SHA、原始组合计数、分片汇总、完成/退出状态和余例映射；主线程重跑这一绑定，并将每片计数与旧完整解压回执逐片再次对照。原始完整回执未因程序更新被冒充重跑。

## 最后证书

positive/cached-recovered-candidate.json 是最后支撑的精确证书。搜索器先由LP选活动支撑，再用有理数解出系数，不能仅凭它宣称成功。主线程实际运行独立的 computation/verify_positive_certificate.py，从支撑重建全部方程，检查185项乘子与30项严格正原子的精确多项式相等。

verify_all_nine34_residuals.py 又实际重核第32轮39个注入证书、6个旧正恒等式及本轮最后1个证书；parent-all46-verification.json 记录全46条及源文件哈希。此前44/36层不属于本包。

## 复现及性能

matrix-builder/build_matrix_cache.py 将每列只转换一次稀疏单项式字典，保持原有基序。两例与旧慢构建逐项相等；最后2381×4958整数矩阵实际3.341秒构建。90MiB可再生数组不归档，元数据和源码保留。旧完整慢构建在产生矩阵前终止，因此没有做该最终矩阵与旧完整路径的SHA比较；独立卷积恒等式验证承担最后认证。

恢复程序只为便于复现增加了输出、缓存和课题路径参数；主线程用参数化版本重新执行并验证结果。不要把数值QR选行或浮点LP当作精确性依据；它们选定后仍需全部整数等式和最终Fraction恒等式检查。

    /Users/zzh/.venvs/main/bin/python matrix-builder/build_matrix_cache.py --I 0,1,2,4,6,9,11,16 --J 0,1,3,4,6,8,11,13,18 --multiplier-degree 2 --atom-degree 2 --cache-dir /Volumes/data/CodexTemp/T2/unfair-01-polynomials/reproduce-nine34-cache

    /Users/zzh/.venvs/main/bin/python positive/recover_from_cached_matrix.py --cache /Volumes/data/CodexTemp/T2/unfair-01-polynomials/reproduce-nine34-cache --topic /Volumes/data/CodexProjects/T2/unfair-01-polynomials --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/reproduce-nine34-positive

    python3 verify_all_nine34_residuals.py --topic /Volumes/data/CodexProjects/T2/unfair-01-polynomials --last-certificate positive/cached-recovered-candidate.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/all46-recheck.json

小收据、完整账本哈希与已执行独立逐项验证相结合才承担有限覆盖；只运行46个余例证书不等于重跑整个搜索。
