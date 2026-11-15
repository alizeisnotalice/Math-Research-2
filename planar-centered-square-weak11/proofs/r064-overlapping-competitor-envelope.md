# R64 supplementary lemma: an overlapping-competitor deletion envelope

2026-11-15. This is a separate general structural lemma; it does not alter the frozen R63 proof or the main R64 mass-polytope proof. It supplies a sufficient condition for a cheap maximum-block deletion, not a proof that all inputs have such a block.

Let A be finitely many distinct positive source atoms for the true centered axis-parallel common-radius square halo at strict threshold1. Let L be the maximum mass of a strictly eligible source subset, and choose a maximum block S, so m(S)=L and its bbox [a,b]×[c,d] has widths w,h<sqrt(L). Define

    L_star(S)=max{m(T): T is a strictly eligible bbox-closed source set,
                        T != S, T intersects S},

with the maximum of an empty collection defined as0. This is a finite source-geometric quantity and satisfies 0≤L_star≤L. The definition includes every such closed eligible set; it does not assert that each one is realized by a positive-area true source state.

## 1. Two actual geometric containers

For a receiver x in H(A), take the union N(x) of the source sets captured by all actual strict witnesses at x. As in the main R64 proof, choosing one witness for each of its finitely many sources and then the largest selected radius produces one actual strict witness capturing exactly N(x). Thus N(x) is itself eligible and bbox-closed.

If x belongs to the true deletion difference H(A)\H(A\S), then N(x) intersects S. Otherwise that actual witness remains unchanged after deletion. There are two cases.

- If N(x)=S, the actual strict witness puts x in the complete-S rectangle

      E_S=(b−sqrt(L)/2,a+sqrt(L)/2)
            ×(d−sqrt(L)/2,c+sqrt(L)/2).

- If N(x)!=S, its mass is at most L_star. Its common strict radius is less than sqrt(L_star)/2, and it captures at least one source of S. Therefore x lies in

      I_star=(a−sqrt(L_star)/2,b+sqrt(L_star)/2)
                ×(c−sqrt(L_star)/2,d+sqrt(L_star)/2).

When L_star=0 the second case is impossible, since every nonempty actual source state has positive mass. It is harmless to retain the displayed possibly degenerate open rectangle as an extra outer set.

Consequently

    H(A)\H(A\S) ⊂ E_S union I_star.                     (1)

The two rectangles have the same center, the center of bbox(S). This is a true-halo deletion containment, not a statement that a nonempty intersection N(x)∩S alone characterizes deletion: residual witnesses may still exist at many of those receivers.

## 2. The universal competitor envelope

The side lengths of E_S are sqrt(L)−w and sqrt(L)−h. The side lengths of I_star are w+sqrt(L_star) and h+sqrt(L_star). Each rectangle is contained in the same centered square of side sqrt(L)+sqrt(L_star). Hence

    Delta_A(S) ≤ (sqrt(L)+sqrt(L_star(S)))².              (2)

There is also a sharper explicitly computable two-rectangle union bound. With u=w/sqrt(L), v=h/sqrt(L), and t=sqrt(L_star/L), (1) gives

    Delta_A(S)/L ≤ (u+t)(v+t)+(1−u)(1−v)
                    −min(u+t,1−u) min(v+t,1−v).         (3)

No additional source masses are paid in these formulas. The old source area is retained in the exact deletion difference on the left. The source-state radius cap merely supplies a smaller outer region for that difference.

**Corollary.** If some maximum block S satisfies L_star(S)≤3L/4, then

    Delta_A(S) ≤ (7/4+sqrt(3))L < (7/2)L.               (4)

Indeed (2) is at most (1+sqrt(3/4))²L=(7/4+sqrt(3))L, and sqrt(3)<7/4 follows exactly from 3<49/16. This argument covers all bbox shapes without a separate thin-axis branch. The empty-competitor case gives the stronger Delta_A(S)≤L.

For a hypothetical counterexample to the R63 maximum-block selection conjecture, **every** maximum block would therefore have to satisfy L_star(S)>3L/4. R63 additionally requires both its coordinate spans to exceed3sqrt(L)/4, because otherwise its original geometric envelope is already below3.5L. These are necessary obstructions to a counterexample, not sufficient conditions for one.

The comparison uses the mass of every overlapping eligible closed competitor as a conservative bound. Computing the exact largest mass among positive-area N-states could lower it further, but would require an additional receiver-state analysis; a closure list by itself does not certify that sharper value.

## 3. Relation to the mass-polytope experiment

The main R64 design constraint bounds all eligible closures other than S by1−eta. It therefore gives L_star≤1−eta, and (2) yields the supplementary design-domain bound Delta_A(S)≤(1+sqrt(1−eta))². For the pilot's tiny eta=1/2048 this remains much larger than its computed deletion cost and does not settle the general target. The continuation code should record L_star as an upper bound from eligible overlapping closures, not label it an exactly reachable receiver-state maximum.

This lemma is analytic. Finite optimization results are not used in its proof. The exact C2 value, a matching legal extremizing sequence, and the assertion that every nonterminal connected input admits a maximum block satisfying the3.5 deletion bound remain open.
