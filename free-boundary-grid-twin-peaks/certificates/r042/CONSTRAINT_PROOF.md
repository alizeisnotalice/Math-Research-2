# Curved-cut intersection DP (finite exact computation)

For a row-width vector h with 1<=h_y<=n-1, put L_h={(x,y):x<h_y}
and R_h=V\L_h. Each is connected because every row meets the fixed left/right
edge. Put d=ceil(n/4). Allowed roots are the vertices of L_h with x<d,
respectively R_h with x>=n-d, that have no neighbor across this cut. E_h is the event
that both induced restrictions have exactly one peak, each in its allowed set.
Every E_h labeling has exactly these two global peaks: a nonroot has a greater
neighbor in its own region, and the roots have no cross-region neighbor. Their
original Manhattan distance is at least n-2ceil(n/4)+1 >= (n-1)/2.
For n>=4 this exceeds n/3, giving a full-n family of linearly separated peaks.
The safe-root condition also forces distance >=3: roots at distance two
would share a neighbor, which would have to belong to both sides. In
particular, every n=5 event here has rho>=3, hence rho>n/2.
The computation is restricted to n<=5; normalization for this infinite family
remains unproved.

The DP reveals vertices in descending label order. For each selected cut and
each candidate next vertex v, let H be its side, A its allowed root set, and
M its neighbor set inside H. If S meets H, require S to meet M. If S does not
meet H, require v in A. The first revealed vertex of H is thus an allowed
root; every subsequent vertex in H has a greater neighbor there. This is
equivalent, in both directions, to E_h. Requiring all selected cuts at every
transition counts their intersection exactly, once per full permutation.

dp[0]=1 and dp[S+v]+=dp[S] for admissible transitions. Masks increase, so a
single increasing-mask pass is valid. All values are at most N! and for n<=5
fit unsigned 128-bit arithmetic. No conditional sampling or mixing assumption
enters. For a family X=sum_h 1_Eh, S1=sum_h #Eh and
S2=sum_h #Eh +2 sum_{h<k} #(Eh intersect Ek) give the rigorous finite bound
P2(X>0)>=S1^2/(A2*S2). Normalized positive mass on an infinite subsequence
is a separate requirement, not obtained from finite computation.

For n=5, h_y in {2,3}, all 32 curves are preserved by simultaneous row
reversal and left-right reflection h->5-h. These graph automorphisms also
preserve the root rules. Thus unordered pairs may be reduced by this group;
diagonal and off-diagonal multiplicities are recorded separately. S1 and S2
must retain these multiplicities. n=3 full 9! enumeration independently
checks all 36 unordered pairs, global two-peak validity, and distance.

For constant h=2 or 3 on n=5 the allowed sets agree exactly with the earlier
straight-cut outer-half root events. Independent old straight-cut DP gives
a cross-check without relying on this program's symmetry reduction.
