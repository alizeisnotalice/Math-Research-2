# Zhang官方源码最小运行

源码/文献核验及证书缺失见references/039-ZeroOnePolynomials-audit.md；本目录只有本次运行脚本及小日志。两个空stdout不是“全范围已证”，而是指定树程序实际运行的结果。

复现：在data临时目录将 https://github.com/dzhang314/ZeroOnePolynomials clone为ZeroOnePolynomials，checkout 582cdcc4a17a84531a5bda6ee5ac8d3fc0ff7766；将repro_minimal.sh置于该clone旁再运行。代码依赖C++20。完整data/proofs未取得、Julia不可用，此脚本不重证≤66。
