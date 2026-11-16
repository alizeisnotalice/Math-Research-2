# MacBook48次GPU搜索与独立重算

Mathematica15.0.0原生GPUArray实际完成24个目标、335544320个数值根组mask，exit0，远端内核16708已退出。40960保存样本全部用NumPy2.5.4独立重构；主线程从归档脚本全量复跑27.722秒，5项生产篡改拒绝。

最大GPU/CPU系数差7.646e-6，最大重乘误差1.284e-5。保存样本中CPU分数<2e-5仅task0/mask131584，为已知二元公平控制；非控制保存样本最小分数0.3391387853。浮点mask全覆盖与保存样本重算不构成精确排除。

完整stdout按2048行压缩成6片，每片<10MiB。先运行 `python3 unpack_stdout.py --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/gpu48-rebuilt.jsonl`，逐片和整流哈希都核验。再 `python3 audit_gpu48.py --input inputs.json --stdout /Volumes/data/CodexTemp/T2/unfair-01-polynomials/gpu48-rebuilt.jsonl --search search.wl --exit exit.json --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/gpu48-recheck.json`。

输入及实际Wolfram stdin程序由SHA绑定；不包括凭据。GitHub官方文件限制于2026-11-16实查 https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github ，本包同时遵守项目更严的单文件10MiB约定。
