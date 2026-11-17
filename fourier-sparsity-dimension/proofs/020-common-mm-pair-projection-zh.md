# 第20轮：同坐标 MM 相位对的两倍翻号块必可投影

状态：主Agent证明，两位独立审查已通过；两机 CPU 与两机 Metal 计算交叉核验完成。一般 A28 仍未解决。

Let X=Y=F2^k, q=2^k, and let
f(x,y)=x·pi(y)+h(y), g(x,y)=x·sigma(y)+j(y),
where pi,sigma are arbitrary permutations. If A=supp(f+g) has size2q,
then A admits an injective linear projection onto a (k+1)-space.
This concerns pairs using the SAME coordinate split; it does not classify
arbitrary bent pairs or arbitrary two-minimum switches.

Put c=# {y:pi(y)!=sigma(y)}, and t=# {y:pi(y)=sigma(y),h(y)!=j(y)}.
Every changed fiber has q/2 points, every full unchanged fiber q points.
Thus c+2t=4. Two permutations cannot differ in exactly one place, and c is
even here, leaving (c,t)=(4,0),(2,1),(0,2).

For (4,0), write S={y1,...,y4}, di=pi(yi)+sigma(yi) !=0.
The multisets of permutation values on S agree, so sum di=0. No odd subset
of the four di sums to0: a singleton is nonzero and a triple sum equals
the remaining nonzero vector. Consequently the linear equations di(v)=1
are consistent (every dependency has even cardinality). Choose such v.
Every four-point subset S of F2^k has a linear quotient to F2² injective
on S: translate one point to0; if the remaining differences span dimension2,
use an isomorphism, and if dimension3 send the three independent differences
to10,01,11. Extend the map to all Y. Its kernel W has dimension k−2 and
W intersects (S+S) only at0. Let K=<(v,0)> direct-sum (0,W).
A difference in K has y-component0, since it is also in S+S; its only
possible nonzero x-component is v, which switches the halfspace membership
in every fiber because di(v)=1. Therefore K∩(A+A)={0} and dimK=k−1.

For (2,1), let a,b be the changed fibers, c the full fiber. The changed
normal is the same nonzero d in both fibers, since the two permutation values
are swapped. Write the two support equations d(x)=ea and d(x)=eb.
Choose lambda in Y* with lambda(c+a)=lambda(c+b)=1; this is possible since
these two nonzero vectors are distinct and hence independent. Set W=kerlambda.
Then delta=a+b is a nonzero element of W. Choose a linear map L:W→X such that
d(L(delta))=1+ea+eb; extend an arbitrary assignment on delta to a basis of W.
Let K={(L(w),w):w in W}. No full/half pair has y-difference in W. A pair
from the two half fibers has y-difference delta and x-difference of d-value
ea+eb, so cannot differ by (L(delta),delta). Within any one fiber, K has no
nonzero vector with y-component0. Hence K avoids all nonzero A differences
and has dimension k−1.

For (0,2), A is the union of two parallel full k-flats, hence an affine
(k+1)-flat. Any complementary (k−1)-space to its direction is a valid kernel.
For k=1 only this last case occurs, and the kernel is zero.

Thus same-coordinate MM pairs satisfy R19's geometric hypothesis in every
dimension, for arbitrary permutations and offset functions. This is wider
than the 293 standard-base eight-variable representatives. It does not
resolve whether every legal two-minimum block has this property.

## 对任意共同扩张的推论

设 f0+V 为全 bent 族，含两两不交块 Ai 的独立指标，每块大小为 q 或2q。
对每个两倍块，只要求存在一个线性坐标分裂，使 f0 与 f0+1_Ai
同时具有上述 MM 形式。不同块允许使用不同分裂。则每个两倍块由本定理
可投影，最小块由第6轮自动为仿射平坦。第19轮主定理遂直接适用：
完整覆盖时 dim V ≤ 9q/2−1，部分覆盖时 dim V ≤ 7q−1+dim V_R。
这里绝不要求 V 的其他方向保持 MM 形式，也不要求各块共用同一分裂。

这一条件仍可能比任意合法翻号对严格。证明并未把一般 bent 相位化为 MM，
也未给一般 A28 的上界或零比值反例。有限跨维验证只能检查上述构造，
不能将这个子类自动提升为全部合法两倍块。
