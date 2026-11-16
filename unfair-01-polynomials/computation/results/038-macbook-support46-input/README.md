# MacBook后继六点45/46层输入与验证原型

搜索用15个Python3.12.13进程，nice10与电脑其他任务共享CPU。源代码通过stdin发送，远端仅RAM，逐I gzip账本流回本机data；无远端大文件。单批14400秒，超时先排空已提交子任务并join，未覆盖即partial_timeout，不算排除。

端到端smoke两chunk共20行精确核验通过。独立verify_support46.py不导入搜索器的分类/规范化函数，逐项重建支持、见证、计数和完整覆盖。最终生产验证器有真实篡改测试；manifest-only永不升级为内容认证，payload逐字重建，生产契约固定45/46、cap6。

复现回归 `python3 test_verify_support46.py`；临时负测写CodexTemp。验证账本 `python3 verify_support46.py <run目录> --report <data上的报告路径> --require-complete`。本轮生产批次仍未完成；parent-live-manifest.json仅检查当时事件与覆盖，没有验证整层内容，严格范围不能升级到46。
