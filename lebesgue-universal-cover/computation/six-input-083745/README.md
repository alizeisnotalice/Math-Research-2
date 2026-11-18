# 0.83745 六输入面积证书

six-input-transfer.md 给出从全域整数核到六个有理点组的证明，outer-083745.md 说明规范域和四层外围。points 目录给出全部 18,028 个点、各组直径及精确面积下端。归档核以 2⁻⁵⁷ 为面积单位，目标整数 120689264254125501。

## 证据范围

生产记录覆盖 1,024 个规范分片和最后的全部 3,273 个闭原叶。普通首次执行 3,071 叶，3,053 叶通过，18 叶达到深度限制。另 196 叶以 .837455 纯五体证书通过。18个深度止叶与另外6个叶均由加密核从整个原叶通过。records/integration-audit.json 是计数、分支、原树地址和转移余量的复核。其 conditional_bound 采用早期 17,513 点，更新的 18,028 点面积下端见 points/normal-gap-083745-audit.json。

四层外围的 8,192 份输出已由独立面积重算生成并逐字节对照。196 个纯五体叶的 348,400 个几何访问已重算，196份回执与196份轨迹均逐字节一致。24 个加密叶的 32,803,666 个节点也已独立几何重算，24份回执的全部数学字段相符。普通核首轮及全体回执已归档，全部普通叶独立重算需要另附其完成回执。本目录的 records/integration-audit.json 只校验既有普通记录。继承的 6,278,533 个原叶排除采用此前各轮的整叶整数分数。native45-prior-index.zip 固定并提供这些分数筛出的全部余叶的位置及未过关的完整整数下端，便于从同一入口重做本靶。历史分数的证明来源见本课题 proofs/tomar-multi-target-scores.md、tomar-cross-index-083743.md 和 rounds/113.md。此前的面积几何运行与本靶的独立重算分别记账。

## 依赖与输入

COMP 取 Math-Research-2 仓库的 lebesgue-universal-cover/computation 目录；本包使用的基线 Git 提交为 d8677d5c4fcd4f6bda8ec64887e1c8747cbac82d。可从其后的提交使用相同的被绑定文件。上游工作区 UPSTREAM 固定为 Tomar 的 10af98606c98ca9ffd0d9f402cae1e2efe74b792，复算时读取该提交的 independent/gen8.py、gen11.py。TREE 所用的 cert11d_08366.tar.xz 压缩包 SHA-256 为 6b5d28bc488edcfde9c6a49f3cb364f9086f445a58f85f166da23727036cf122，解压到以 SHA256SUMS 和 part0000.json 等文件为直接子项的目录。原树只提供分割坐标。 [Tomar 补充材料](https://doi.org/10.5281/zenodo.23717912) · [冻结上游源码](https://github.com/tejasvi/lebesgue-cover-bound/tree/10af98606c98ca9ffd0d9f402cae1e2efe74b792)

复算环境为 Python 3.13.18、NumPy 2.5.5、python-flint 0.9.0、SciPy 1.18.1 与 Apple Clang 21.0.0。面积整数核使用 C++17、带符号溢出陷阱；几何区间采用160位 Arb。使用启用断言的 Python。构造程序使用确定参数，无随机种子。数值提案只在确定性有向核验通过后计入证明。

## 复算

以下以包目录为当前目录。COMP、TREE、UPSTREAM、WORK、CACHE 均由调用者指定实际路径；WORK 为新建工作目录，CACHE 为可写 Clang 模块缓存。首先核对全包文件散列。macOS 可用 shasum -a 256 -c SHA256SUMS，GNU 系统可用 sha256sum -c SHA256SUMS。各分卷另有逐项散列清单。

```sh
python -I -B prepare_native45.py "$COMP" "$TREE" "$WORK/prepared"
T32="$WORK/prepared/table32.npz"
T64="$WORK/prepared/table64.npz"
S="$WORK/prepared/seeds"
python -I -B audit_outer_083745.py seed "$COMP" "$TREE" "$S" "$WORK/prepared/index" "$T32" "$T64" volumes "$WORK/screen" > "$WORK/seed-audit.json"
python -I -B "$COMP/replay_tomar_signed_translation_screen.py" "$COMP" "$TREE" "$T32" "$S" "$WORK/prepared/seed.json" "$WORK/screen" "$CACHE" --workers 4
python -I -B audit_outer_083745.py screen "$COMP" "$TREE" "$S" "$WORK/prepared/index" "$T32" "$T64" volumes "$WORK/screen" > "$WORK/screen-audit.json"
python -I -B integrate_native45.py volumes "$COMP" "$TREE" "$WORK/integration.json"
python -I -B pure45.py volumes "$UPSTREAM" "$TREE" "$CACHE" "$WORK/pure" "$WORK/pure.json" --workers 4
python -I -B replay_dense_083745.py "$COMP" "$TREE" "$T64" "$WORK/dense" "$CACHE" "$WORK/dense.json" --workers 4
python -I -B replay_ordinary_083745.py "$COMP" "$TREE" "$T32" "$T64" "$WORK/ordinary" "$CACHE" --workers 4
python -I -B transfer45.py "$COMP" "$TREE" "$T32" "$T64" volumes/native45-receipts.tar.gz "$WORK/transfer.json"
python -I -B certify_points_083745.py "$COMP" "$TREE" "$WORK/points"
```

独立复算时为各面积核指定全新的输出目录。加密核允许用 --check-existing 仅核对先前已实际执行的结果；此模式不产生新的几何计算。pure45.py 支持接续已落盘的叶。比较数学回执时，编译二进制的 UUID 与软件版本元数据须与逐叶数学内容分别处理。六组点表 SHA-256 为 f78e61ab6e365236e0962cdfac9c040d2a0c61b2cf794be5fada36860e39ccfc。
