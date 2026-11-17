# R036 dyadic-box 核验运行收据

- 日期：2026-11-17（Asia/Shanghai）。
- 解释器：`Python 3.15.0`。
- 命令：

  ```sh
  time nice -n 10 python3 \
    verify_dyadic_shell_samples.py
  ```

  退出码为 0；实际计时为 `0.31s user 0.01s system 64% cpu 0.495 total`。
  输出保存为 `dyadic-shell-sample-check.json`。其中固定种子下，$n=16$ 实际检查
  dyadic 块 $m=2,4$，$n=32$ 实际检查 $m=2,4,8$；每条样本均以 Python `assert`
  核对壳层、唯一早邻点、回插秩界和块和界。

- 整数恒等式单独以

  ```sh
  nice -n 10 python3 -c 'import math; assert all(
      max(0, math.floor(math.log2(n//2))-1)==max(0, n.bit_length()-1-2)
      for n in range(2,1_000_001)); print("pass")'
  ```

  核对 $2\le n\le10^6$，输出 `pass`，退出码为 0。该有限测试只复核实现；文中
  恒等式本身由二进制位长直接证明。

- 运行结束后未留下上述 Python 进程。没有使用 SSH、GPU 或大内存计算。

冻结对象的 SHA-256：

```
be382526fb0886d846dfebc47b0d354b96c58c01d5821e250504d9847fdb65a4  DYADIC_BOX_PROOF.md
7130a13dd58bdf377de8228b6e6a92e4fc695453a1e63705bac25d235b0d6484  verify_dyadic_shell_samples.py
a496b5414b37129188e20716e55dd6d6eccdacb73080ca15197bc262406e70df  dyadic-shell-sample-check.json
```
