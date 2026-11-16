# H8(x³)H8(x⁵) 精确复现

Python3.13.18，验证只用标准库整数/Fraction。生成中心另用mpmath110位。证明见 `../../../proofs/029-H835-exact-theorem.md` 和独立解析证明；本包不解决任意二元多项式。

从本目录执行，先确认data挂载；只在指定临时目录解压，不改归档证书：

```sh
export H835_INPUT_DIR=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/H835-reproduction
python3 - <<'PY'
import gzip, os
from pathlib import Path
assert Path('/Volumes/data').is_mount()
p=Path(os.environ['H835_INPUT_DIR']); p.mkdir(parents=True, exist_ok=True)
for name in ['H8-a3b5-canonical-masks.u32le','witnesses.u8']:
    (p/name).write_bytes(gzip.decompress(Path(name+'.gz').read_bytes()))
PY
python3 verify_h835_certificate.py --balls H8-a3b5-root-balls.json \
  --ledger "$H835_INPUT_DIR/H8-a3b5-canonical-masks.u32le" \
  --certificate certificate.json --witnesses "$H835_INPUT_DIR/witnesses.u8" \
  --galois H8-mother-certificate.json --record "$H835_INPUT_DIR/verification.json"
python3 test_h835_u8_guards.py
python3 test_frontend_semantic_tampers.py
python3 test_reverse_ball_kernel.py
```

最终主线程完整回执为 `parent-ball-verification.json`，绑定最终验证器、前端、倒序核、Rouché helper及模Galois核的SHA。`agent-pre-refactor-verification.json`仅保留旧版本实际全量记录，不拿它替代最终版本认证。

严格验证不信任搜索状态。`certificate.json`绑定根盘和未压缩u32le流；每个二字节见证只指定因子侧和距最高次的偏移。验证器自己重建规范流和系数区间。标准项必须精确对应117901063；其余每项严格负才接受。

搜索重跑须用新空目录：

```sh
python3 search_top.py --balls H8-a3b5-root-balls.json \
  --ledger "$H835_INPUT_DIR/H8-a3b5-canonical-masks.u32le" \
  --out "$H835_INPUT_DIR/new-search"
```

完整根盘生成器和独立前端一并保存。生成器写当前目录的root-balls文件，应复制到新的临时目录执行。测试只作了本包路径/`H835_INPUT_DIR`适配；数学验证代码保持冻结副本。
