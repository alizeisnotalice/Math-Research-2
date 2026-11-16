# r024 树计数冻结回执

冻结程序：[tree_peak_dp_r024_frozen.py](tree_peak_dp_r024_frozen.py)，SHA-256
`76eddd4aa26bfe45e8e736adf95413317197f71db780534aad639050f5ff529b`。

`rerun_frozen.py` 实际自然结束，输出
[frozen-rerun.json](frozen-rerun.json) SHA-256
`f967ef59c53ac397c72094afb4015c435f546e5602a42c4f7844e085b51b8d49`。
该次重算覆盖单满二叉树 \(h=1,\ldots,9\) 和两棵等高满二叉树根—根连接的树
\(h=1,\ldots,7\)，均为精确有理数。资源回执：墙钟429.30秒、用户CPU252.03秒、
最大RSS 25,935,872字节、交换0。

独立 hook/cut 验证器 [verify_cut_formula.py](verify_cut_formula.py) 的最终 SHA-256
为 `88b6d164dd33d9b80031d8ab7d16b4a950a296268e3364e9ac95f7578079f718`。实际检查
P4、P5、满二叉h3、双块h2/h3/h4：每例的切边公式整数 \(B\) 与积分消息输出相等；
并对每条相对于指定根的切边检查含祖先 hook 修正因子 \(G\) 的分解。原先省略 \(G\)
的分解式已经撤回，见 [CUT_FORMULA.md](CUT_FORMULA.md)。

有限结构反例：双块h3（30点、最大度3）满足
\(A_2/B=838907595/253657849>3\)，严格否定候选 \(A_2\le\Delta B\)。这是一个
辅助 normalizer 不等式的反例，不是自由方格 twin-peaks 原题的反例。双块h1…h7
比例上升至约17.5825，但尚无全高度下界证明，因此不写成 \(A_2/B\) 发散结论。
