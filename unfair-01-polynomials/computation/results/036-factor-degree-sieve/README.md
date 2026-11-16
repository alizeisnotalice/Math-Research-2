# 3190个既有未知有理因子的证书补齐

固定源为 ../033-qfactor-certificates/factor-certificates-d20-d24.json.gz，SHA c0527302d8d9612b1b4397a29ef038dd5718a5b2a3bdaa326b520cc906466410。
9进程发现器已结束；3183个次数筛证书、3个整余式证书，4个仍未知。输入启动快照 input.json 的 running 是历史值，完成状态见 result.json。
独立整数验证不调用SymPy；完整源绑定和六类语义负测见 bound-verified.json。数学依据见[证明](../../../proofs/036-factor-certificate-completion.md)。

从本目录复验：

```sh
mkdir -p /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay036/factors
gzip -dc rows.jsonl.gz > /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay036/factors/rows.jsonl
/Users/zzh/.venvs/main/bin/python verify_bound_and_lifts.py --source ../033-qfactor-certificates/factor-certificates-d20-d24.json.gz --rows /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay036/factors/rows.jsonl --lifts lift-eisenstein.json --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/replay036/factors/audit.json
```

Python3.13.18、SymPy1.14.0仅用于发现。独立验证7.255秒通过；不会因此获得全部实根分配排除。
