# Whole-leaf scores at 0.83743

The reduction uses the same 2,048 original roots and the same canonical slice of 6,710,314 closed leaves as the 0.83742 proof. Scores are whole-leaf lower bounds in units of $2^{-57}$. A score passes a rational target $t$ when the integer score strictly exceeds $2^{57}t$. Discarded scores in each input index reach the truncation level 0.83745. For every queried target at most 0.83745, intersect the two lists of retained original positions and take the larger of the certified scores at each common position. The tree manifest, canonical slice, root and radial certificate are checked for identity before this operation.

The first index is the signed four-layer index at 0.83742. The second uses the independent standard four-layer run at 0.83741 and the published complete index at 0.83740. Replaying its packet chain from the seed reproduces all 4,096 packet hashes. A final maximum uses the 2,335 exact integer witness records from the previous GPU-guided search. Only their strict CPU-verified integer scores enter the comparison.

| Target | Signed index | Standard index | Both indices | With integer witnesses |
| --- | ---: | ---: | ---: | ---: |
| 0.83743 | 268,461 | 420,587 | 255,973 | 255,973 |
| 0.83744 | 431,780 | 517,568 | 419,095 | 419,091 |

Starting from the two authenticated index directories, these commands regenerate the input subsets and the four signed stages. The optional argument `--replay-standard CERT TABLE32 CLANG_CACHE` on the first command recomputes the standard geometric trees and compares their complete packet outputs.

```sh
python -I -B computation/replay_tomar_cross_index_083743.py computation INDEX_083740 INDEX_083742 WORK/index
python -I -B computation/replay_tomar_signed_translation_screen.py computation CERT TABLE32 WORK/index/seed43 WORK/index/seed43.json WORK/signed CLANG_CACHE --workers 4
python -B computation/audit_tomar_083731_score_layers.py WORK/signed WORK/index/seed43 WORK/index/seed43.json WORK/index/seed44 WORK/index/seed44.json --next-emit WORK/next
```

The signed kernel certifies the original closed box by combining its clipped score with the radial lower bound for the removed part. The portable input reconstruction checks the common tree identity, the saved four-layer packet chain, the authentication of both input indices, and the exact output hash of each new seed. Each surviving original leaf still requires a terminal proof.
