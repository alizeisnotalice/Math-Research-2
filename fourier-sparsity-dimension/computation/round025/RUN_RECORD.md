# 第25轮实际验证

本机Python3.15.0、Apple M6、macOS27.2、Swift6.4.2；MacBook Apple M5 Pro。所有数学检查采用整数/F₂运算，无浮点容差。

sidon_seed_audit.py构造k2..6的有限域，逐无序对验证评价和非零且互异，并独立检查立方差分每纤维大小2。k2完整65536函数枚举得到896个bent，384个兼容H的函数，逐纤维均为指定split的MM。

两机Swift各直接计算65536×16=1048576个Walsh系数，CPU用16×16定义式字符和，GPU用四阶段整数FWHT，零分歧。另以纤维识别器及24个置换×16个偏移的独立构造核对384个MM函数；两机完整384元素集合与Python完全相同，不仅比较计数。

weight_four_obstruction.py对k2,3,4分别完整检查8、128、32768个g+H相位，全bent；指定w=0纤维的z谱大小4，确非固定split MM。该计算验证四点奇偶缺失的具体负例；全维充要判据由解析证明承担。

中间Swift源的Metal内核buffer属性漏闭合括号，导致运行期MSL编译失败；主Agent修正后本机完整执行通过，远端亦采用修正版重新运行。失败输出不作为最终证据。最终源及输出绑定、编译/运行命令和退出码见各receipt；不将一次SSH观察超时解释为仍运行或已终止。

复验存档运行`python3 verify_receipts.py`。Python脚本会覆写同目录JSON，建议在副本复现；Swift示例：

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-module-cache R25SidonMetal.swift -o /Volumes/data/CodexTemp/T2/a28/r025-local/reproduce
/Volumes/data/CodexTemp/T2/a28/r025-local/reproduce /Volumes/data/CodexTemp/T2/a28/r025-local/reproduce.json
```

没有打开MATLAB，无持续满载声明。种子和字段的有限核验不替代一般A28证明。
