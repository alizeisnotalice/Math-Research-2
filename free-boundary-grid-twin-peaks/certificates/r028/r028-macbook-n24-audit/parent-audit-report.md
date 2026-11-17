# R028：MacBook 返回的 n=24 完整批次独立审计

## 范围和证据等级

审计对象为只读原件
`raw/n24-p1048576-r32-l16.jsonl.gz` 解压后的旧式
uniform-move、$n=24$、$P=1048576$、$R=32$、$L=16$、域
`9f29010000000000` 的单一有限批次。本审计没有 SSH、没有重跑 GPU、没有修改
原件。它验证保存的终端原始输出中的整数恒等式、标识、冻结聚合可复现性和存储小数上的
必要条件插件检查；它**不**证明 proposal 已混合、也不把任何比例或 SE 当作
原条件律的置信区间。

这里的生产批审计不是全 rank 轨迹的逐步独立重放：证据对象仅为保存的终端
group 记录及其冻结聚合。

## 完整性与运行记录

原始 JSONL 的 SHA-256 为
`c0695ddbaf7d43494bc0169323cf04c795dfff41aa960cea71bbf5228161ccc5`；
Swift 源的 SHA-256 为
`01b77bca1735f07293eae5ffeea0326810fdc5ea785c7a9e99d9030cba9b0702`；它们都与
sidecar 相符。二进制哈希为
`137aa827fd2fcfbacd9a40aaa7059ee9dd8f38a3616e0c849724d2143f10f2cd`。

对保存的终端原始 JSONL 作流式审计验证：恰有 1 条 metadata、32 条 replicate（ID 为 $0,\ldots,31$）、
1 条 run_end；每个 seed 精确等于

$$
(\mathtt{0xA022000000000000}+\mathrm{domain}+24\cdot10^6+2\cdot10^4+r)
\bmod2^{64}.
$$

所有 32 个 seed 唯一，metadata 与 sidecar 的参数完全一致，run_end 的
`wall_seconds` 有限。driver 与 trace exit 均为 0；driver 记录的总 wall time 为
46796.18367385864 s，run_end 为 46789.3505589962 s。资源记录显示峰值 RSS
9021259776 B、峰值 footprint 11970108808 B、零 swap。`driver.stdout` 和
`driver.stderr` 为空；本审计未自行推断远端 PID 状态。

## 整数载荷与距离直方图

每条 replicate 有 4096 个 256-particle group，共 $32\times4096=131072$ 组。
逐组验证了

$$ LF-2LG=4LD $$

（记录中的残差字段均为零），以及每个距离行的整数和严格等于该组的 $L$。
所有距离行长度为 $2n-1=47$，第 0、1 bin 为零。该终端原始输出所有保存 group 的整数合计为

$$\sum L=69913412,\qquad \sum L^2=228067664,$$

且独立累计的 47-bin 全局直方图和恰为 69913412。该结果证明已保存的整数
payload 自洽；它不验证 SMC 权重的统计代表性。

## 冻结聚合器重放

使用项目冻结的
`r024-logscaled/frozen_r023_aggregate_tau2_distance.py`、sidecar 和原 Swift
源，分别以 $\varepsilon=1,1/2,1/4$ 本机重放。三个 replay JSON 与原件的
`aggregate-epsilon-*.json` 都逐键结构相等（不是近似比较）。

| $\varepsilon$ | replay SHA-256 |
|---:|:--|
| $1$ | `ce82dd69de45b3075e3737b4a19505e09e47c0686dffdd959867b9e95c70ca09` |
| $1/2$ | `2147e4ebed52c8cf01594d40fe46c8aa1e951c731739bf37fb49d24a275a92dd` |
| $1/4$ | `08caf7c44816a3ac61d709227cb33ccdc969a6787db225dc5662e7b5d95667d7` |

重放的有限插件量为

$$\widehat{B/A}=0.8592774459438061,\qquad
\widehat Q=0.2785758124496006.$$

严格距离 tail（聚合器以 `Fraction` 交叉乘法确定阈值 bin）为：
$\widehat{\mathbb P}(d>24)=1.1281701\cdot10^{-11}$，
$\widehat{\mathbb P}(d>12)=0.99999965218$，
$\widehat{\mathbb P}(d>6)$ 在存储精度下为 1。这些仍是该高度集中加权批次的
有限插件值。

## R025 必要 cone：精确符号结果

对重放的 $\varepsilon=1/2$ JSON，所有十进制字符串按 `Fraction` 精确读取并先
归一化直方图。R025 线性必要条件的两边为

$$
\widehat{B/A}=0.8592774459438061,\qquad
\sum_{d\ge2}\frac{d}{2(d-1)}\widehat h_d=0.5299228948540377.
$$

因此精确有理数 slack 是负的，数值为约 $-0.3293545511$；本批次的线性 cone
插件**违反**该必要条件。阈值 tail 插件在 $r=2,6,12$ 也违反，$r=18$ 通过。
精确分子/分母、每个阈值的符号均保存在 `raw-and-cone-audit.json`，因此这个
“负号”不来自浮点舍入。

这不是原题的反例，也不能单独宣布程序错误：它是一个高度集中重要性权重批次的
失败有限诊断。它足以阻止把本批次当作 $B/A$、距离律或远峰概率的可信数值证据。

## 确定性二择误差认证

令保存的距离直方图十进制小数以 `Fraction` 精确读入并重新归一化为
$\widehat q$；令 $g(d)=d/[2(d-1)]$（$2\le d\le46$）。则

$$\operatorname{osc}(g)=1-\frac{23}{45}=\frac{22}{45}.$$

前节的精确锥差为
$\delta=\widehat{B/A}-\mathbb E_{\widehat q}g>0$。若真实对象满足同一条
必要 cone，则三角不等式和（这里 $\operatorname{TV}$ 采用标准
$\frac12\ell^1$ 约定）
$|\mathbb E_{\widehat q}g-\mathbb E_qg|\le\operatorname{osc}(g)\,
\operatorname{TV}(\widehat q,q)$ 给出完全确定性的结论

$$
\max\left\{\left|\widehat{B/A}-B/A\right|,
\operatorname{TV}(\widehat q,q)\right\}
\ge\frac{\delta}{1+22/45}=\frac{45\delta}{67}
\approx0.22120828058268027>\frac{11}{50}.
$$

比较 $11/50$ 亦以精确有理数完成。此认证不使用置信区间、独立 replicate
近似、ESS 模型或混合假设，且不归因于任何特定 bug；它只说明在必要 cone 成立
的前提下，保存插件的 $B/A$ 与保存距离律至少有一者不可能同时接近真实对象。
详细分子分母见 `cone-discrepancy-certificate.json`。

## ESS 限制

冻结聚合器报告的 replicate-weight 诊断为

| 权重 | ESS | 最大权重份额 |
|:--|--:|--:|
| tau=2 normalizer | 1.0000151157654955 | 0.9999924421923426 |
| lift normalizer | 1.0000391906436639 | 0.9999804051239098 |

这说明几乎全部估计质量落在单一 replicate；输出的 delta-SE 不能解释为通常的
独立 replicate 误差，亦没有混合或准确度保证。上述 cone 失败与这种集中相容，
但两者之间没有由本审计证明的因果结论。

## 复现

```sh
python3 -B frozen_r023_aggregate_tau2_distance.py RAW --sidecar SIDECAR \
  --source tau2_distance_terminal_r023.swift --epsilon 1/2 > replay-epsilon-1_2.json
python3 -B audit_raw_and_cone.py RAW SIDECAR tau2_distance_terminal_r023.swift \
  replay-epsilon-1_2.json raw-and-cone-audit.json
```

`audit_raw_and_cone.py` hash：
`3ea781e502bc4a08734cd8625ee2fa5147a99369bc3699c5ef26ed29c1b08b62`；
其输出 hash：
`f2ea8b5ead3374640200052aa8d56a28792843d48d412fec12735389d37da319`。
确定性误差证书脚本 hash：
`dc92b4bd819d968d4d8cf7648c3404115acaa06d389615f9372c700c40cb4e00`；
输出 hash：
`398c7e175f7b65766877af2b4f15fcb131c52eef5b20ab527e9e92e3806921dc`。
