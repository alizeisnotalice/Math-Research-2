# H7(3,4) 固定目标完整精确证书

12,158个必要根布局中12,157个有严格负系数，余下mask6407047精确为h(x³)/h(x⁴)。这不是一般猜想的证明。搜索器和验证器独立：前者128bit有向dyadic全系数，后者从根盘重新构造矩形盒并用Fraction截断到指定系数。

从本目录运行：
```sh
python3 check_h7_target.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H7-target.json
python3 verify_h7_certificate.py --balls H7-a3b4-mixed-root-balls.json --ledger H7-a3b4-mixed-layout-ledger.json --certificate H7-certificate.json --galois H7-mother-certificate.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H7-recheck.json
python3 test_frontend.py --balls H7-a3b4-mixed-root-balls.json --ledger H7-a3b4-mixed-layout-ledger.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H7-front-tests.json
python3 test_certificate.py --balls H7-a3b4-mixed-root-balls.json --ledger H7-a3b4-mixed-layout-ledger.json --certificate H7-certificate.json --galois H7-mother-certificate.json --scratch /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H7-cert-tests
python3 generation/search_h7_dyadic.py --balls H7-a3b4-mixed-root-balls.json --parent-masks canonical-masks.json --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H7-new-search
```

搜索输出路径要求为空；每100布局和finally原子保存检查点，异常时保留partial状态。搜索器只把canonical-masks.json作独立集合对照，仍自己生成所有选择。重新生成根盘可在独立临时目录运行generation/generate_h7_mixed_balls.py（需mpmath）；不必重新生成数值中心才能检查现有精确根盘。

实际Python3.13.18；验证仅标准库。full-verification.json为主线程完整Fraction结果；negative-tests.json为11项前端回归，coefficient-negative-tests.json为8项系数语义回归。review-*是另一次前端/S7复跑。搜索lo/hi是建议界；检查器不凭它们判符号，而用自身Fraction系数上端给出结论。96bit余量只是将独立严格上端向上舍入后的保守正下界。

完整49个互不相交目标根圆盘，同时证明两族平方自由互素。生成器/搜索器同名helper位于generation，验证使用自己的h7_exact_common.py；两者不混用。SHA256SUMS覆盖本包最终文件。不存在未经认证的端点pin或数值成功标志替代验证。
