# R19：重量恰为两倍最小值时的文献边界

检索日期：2026-11-17。目标是：对 \(RM(r,n)\) 中重量 \(2^{n-r+1}\) 的支撑 \(A\)，是否必有 \((r-1)\)-维线性核 \(K\) 满足 \(K\cap(A+A)=\{0\}\)。以下资料不能代替该命题。本轮解析预算不依赖此类外部归约。

- Borissov、Manev、Nikova，*On the Non-Minimal codewords of weight 2dmin in the binary Reed-Muller Code*，WCC 2001，103–110。[作者机构公开稿](https://www.esat.kuleuven.be/cosic/publications/article-413.pdf)，PDF 第2–5页。其两不交最小仿射平坦分解要求码字**非极小**；没有覆盖任意重量2d码字。PDF SHA-256：`12310bb1996d6c02abba1cd5ad6aac4652c918fb1f9f275875968fa2e705bff1`。
- Klein、Storme，*On the non-minimality of the largest weight codewords in the binary Reed-Muller codes*，AMC 5 (2011), 333–337，[作者库全文](https://backoffice.biblio.ugent.be/download/2140143/2939938)，DOI 10.3934/amc.2011.5.333。第3页 Lemma 1 在 \(r\ge4\)、重量严格小于 \(2.5d\) 时给出两个独立超平面覆盖，包含重量2d，但未给所需避差线性核。PDF SHA-256：`5731b8b717d8aa208ae7a5c79be1524018cb193c2eb8322f0895361b76b7c569`。
- Schillewaert、Storme、Thas，*Minimal codewords in Reed–Muller codes*，DCC 54 (2010), 273–286，[2009作者扩展摘要](https://cage.ugent.be/geometry/Files/342/MinimalRevision.pdf)，DOI 10.1007/s10623-009-9323-x。第2页 Theorem 1.8 的正规形要求重量**严格小于2d**，不能用于等号。PDF SHA-256：`4daeb3b65e41816e780532388489f65dd7a861952b35b570bd1aee6e34369978`。

上述 Kasami–Tokura–Azumi 归属来自已读论文中的陈述；未获得其1976原文可读全文，不能称已核验原证明。子Agent读取上述公开稿，主Agent复核相关段落及机构链接；全文仅存临时目录，不纳入公开仓库。本轮未找到目标一般可投影定理或反例，不把有限检索阴性当作不存在文献的证明。
