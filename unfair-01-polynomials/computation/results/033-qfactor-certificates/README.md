# 有理因子清单的部分独立不可约性认证

本包绑定第32轮 search-d20-d24 中16,339个不同的有理因子。独立纯整数验证器实际重算通过13,149个：一次因子1个、分圆多项式32个、素数模约化Rabin证书13,116个。另3,190个保持 unknown；其中3,251个母式、9,716个分配涉及 unknown，因此不能声称整个有理清单已经完备认证，更不能据此排除任意实因子。

原始证书16,528,417字节以确定性gzip归档；原始与压缩SHA均保存在 manifest.json。这些是本轮生成的精确整数证书，不是从历史工具输出推定。

主审发现早期验证器没有检查模数素性、次数绑定、Rabin检验项完整性及重复证书行，已撤回早期验收回执。这里只归档加固后的最终程序与实际重跑回执（162.811秒）。它独立试除判素、重算次数所有不同素因子、要求完整无重复的gcd列表，并拒绝不规范系数与重复行。生产负测实际拒绝删gcd、复合模数、错误次数和重复行。生成和复核进程均已退出。

判据依据：对首一次数 \(n\) 的有限域多项式，\(x^{p^n}=x\pmod f\) 迫其不同不可约因子的次数都整除 \(n\)；对每个素数 \(r\mid n\)，再要求 \(\gcd(f,x^{p^{n/r}}-x)=1\)，便排除每个真因子次数。因此 \(f\) 不可约。首一整数多项式的同次模素数不可约约化由Gauss引理推出有理不可约。分圆项以整数递归复算 \(\Phi_n\)，并使用分圆多项式在有理数域不可约的标准定理；不把仅相等于该递归表达式说成重新证明该标准定理。

复现时先用第32轮的 unpack_ledgers.py 在 data 临时目录解压 h-coverage.jsonl。将本包gzip解压到同一临时工作区，并核对 manifest.json 中原始SHA，然后运行：

    python3 scripts/verify_factor_certificates_integer.py --search-out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round032-reproduce/search-d20-d24 --certificates /Volumes/data/CodexTemp/T2/unfair-01-polynomials/factor-certificates-d20-d24.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/factor-certificates-recheck.json

验证器不调用SymPy分解。生成器使用既有SymPy环境；搜索分配和负系数见证的完整性仍需第32轮分配验证器承担。本包不可约性检查只针对绑定清单的各个因子，不独自复证所有母式与分配。
