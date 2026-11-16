# R77: a qualified uniform-grid obstruction to the total-local-mass mixed proxy

2026-11-16. The maximal basis consists of centered open axis-parallel squares with one common radius and strict threshold1. The source measures are finite positive atomic. This result disproves a proposed auxiliary payment, not the actual uniform-grid deletion bound or a bound for C2.

## 1. The qualified input and the precise proxy being tested

Let n=20, e=1/400, h=21/400, w=399/400, and let sigma place mass e at every point of S={0,h,...,19h}². Put

    a=71/200, d=171/200=a+1/2,
    q1=(-a,0), q2=(0,-d),
    nu=(1/4)delta_q1+(1/2)delta_q2,
    Gamma=[-1/2,-e]², P=(-1,1/2-e)².

There are100 internal atoms in P, with total mass m=1/4. Both outside atoms are in P. The original402-source input satisfies every closed-square cap

    (sigma+nu)(Q_R)<=max(1,R²),                     (1)

and S is its unique maximum eligible block of mass1. It is connected and nonterminal in the R63 sense. These qualifications are established in Section2 without replacing internal positions.

Define alpha=m delta_0 solely for constructing the following auxiliary receivers. For i=1,2, let A_i be the set of centers for which the label consisting of alpha and the i-th outside atom has a strict witness with r<1/2. Let A_12 be the analogous receiver of alpha and both outside atoms. The proposed mixed-only proxy is

    M=[(A_1 union A_2 union A_12)\H(nu)] intersect Gamma. (2)

The outside masses and positions in every label are the same actual nu, and H(nu) is its complete residual halo over all radii. The explicit alpha-singleton receiver is not added to the union. Points lying in that singleton receiver and also in a mixed receiver remain in(2); the singleton receiver is not subtracted either.

We prove

    |M|=43sqrt(3)/200-101sqrt(2)/1600-24641/160000
       >206483/1600000>1/8=m/2.                    (3)

The same lower bound also exceeds the stronger single-corner sufficient threshold1/8+3/(4n⁴)=80003/640000, with strict gap12951/3200000. This only rules out using this proxy for that single-corner payment; it does not imply failure of the four-corner sum, since other corners may have slack and independently pasted inputs need not remain feasible.

Thus original-input all-scale feasibility, the actual corner/source windows, the radius cutoff and full residual subtraction do not suffice to make this particular total-local-mass proxy obey half-mass payment. The original positions of the captured internal prefixes still matter.

## 2. Original full-grid qualification

The unequal-axis cap lemma in [the companion proof](r077-mixed-only-obstruction-and-unequal-axis-caps.md), Section3, gives an exact finite qualification. At these effective positions, one must check n+n+n²=440 rational rows:

    b1+k/n <=max(1,[a+(k-1)h]²),
    b2+l/n <=max(1,[d+(l-1)h]²),
    b1+b2+kl/n² <=max(1,max[a+(k-1)h,d+(l-1)h]²).

The archived production certificate `computation/results/round077-mixed-proxy/n20-qualification.json` evaluates all440 rows for the stated rational input. Their minimum slack is16000/160000=1/10. By the companion lemma this proves(1), not merely unit-square feasibility.

The same production run separately reconstructs all45302 nonempty tight source-bbox closures of the actual402-source input. Every cap passes;7386 closures are strictly eligible. The maximum eligible mass is1, attained only by S. This is a complete finite closure check, not enumeration of all subsets of402 sources and not a halo integration. The strict prefix-row slack also explains uniqueness: every label containing an outside source has strict cap slack, so it cannot be both eligible and have mass at least1; a proper old-only subset has mass below1.

Connectivity has a short independent certificate. All grid neighbors form eligible two-point edges since2e>h². The origin and q1 form an eligible pair since e+1/4>a². The outside pair connects q2 because3/4>d². These edges connect all402 sources. The whole-source span is w+d=741/400, and

    (w+d)²-(35/33)(7/4)=8319673/5280000>0.          (4)

Hence the input is nonterminal. The activity menu of all100 original old prefixes with both outside atoms contains9 strict active rows and91 null rows. The active rows are exactly(k,1), k=1,...,9, so their largest captured internal mass is9/400, far below the collapsed100/400. This is a statement about strict joint full labels, not their surviving deletion area. This places the example outside the R76 joint-null hypothesis; that count alone makes no assertion about actual corner deletion area.

## 3. Complete mixed and residual receiver geometry

Use X=-x,Y=-y, so the corner becomes Gamma+=[e,1/2]². Set

    u=sqrt(2)/4, v=sqrt(3)/4.

The elementary orderings used below are

    1/4<a-e<u<a<v<1/2,
    d=a+1/2<2v, d-u>1/2.                           (5)

All receiver formulas below are stated with closed endpoints for area; strict witnesses omit their appropriate boundaries, which have zero area. The mixed radii are u,v,1/2, respectively. The last value is a supremum: actual witnesses always have r<1/2. Direct intersection of the containing-center intervals gives

    A_1 intersect Gamma+ = [e,u]²,
    A_2 intersect Gamma+ = [e,v] x [d-v,v],
    A_12 intersect Gamma+ = T=[e,1/2] x [a,1/2].     (6)

The complete residual has three labels, the two singletons and their pair. Their radius thresholds are1/4,u,v. The second singleton misses Gamma+ by d-u>1/2. The residual pair has precisely the same corner receiver as A_2, because m=b1 and the first outside constraint is already inactive at radius v on this crop. Consequently A_2 contributes nothing to M. The first singleton intersects only A_1, in a rectangle of area

    (u-a+1/4)(1/4-e).                               (7)

The pair residual lies inside T, with area

    (v-e)(2v-d).                                    (8)

It does not meet A_1. Indeed d-v>a>u, whereas the first singleton has Y<=1/4<a and so does not meet T. The two surviving full-receiver pieces A_1 and T are disjoint since u<a. Thus no overlap or residual label is missing from

    |M|=(u-e)²-(u-a+1/4)(1/4-e)
          +(1/2-e)(1/2-a)-(v-e)(2v-d).              (9)

Expanding(9) gives(3). This is an exact area calculation for the auxiliary proxy, with the same nu throughout; it is not a sum of independently chosen one-atom optimizers.

For an elementary strict comparison, use sqrt(3)>433/250 and sqrt(2)<283/200. Both follow by squaring positive rationals. The coefficient of sqrt(3) is positive and that of sqrt(2) is negative, so(3) exceeds206483/1600000. The excess over1/8 is at least the strict rational amount6483/1600000. The displayed decimal0.129112442502507 is only a convenient approximation of the exact radical expression.

## 4. A positive-area separation from the actual full halo

The discrepancy is visible without integrating the actual402-source halo. Let

    E=[-49/100,-47/100]², |E|=1/2500.               (10)

For every x in E, the collapsed origin and both outside positions lie in the open square with common radius99/200. Their collapsed total mass is1, and

    1-4(99/200)²=199/10000>0.

Thus E is contained in the proxy triple receiver, with the required radius strictly below1/2. Every outside atom has one coordinate equal to zero, so a square centered in E capturing either atom requires r>47/100. All residual mass is3/4<4(47/100)². Hence E avoids H(nu), and E is contained in M. It also lies outside the collapsed old singleton halo of radius1/4.

For the actual full input, radii r<=47/100 capture no source. If47/100<r<1/2, the upper endpoints of the receiving square are less than3/100<h. At most the actual origin can be captured from S, along with the two outside atoms. Their total mass is at most

    e+3/4=301/400<2209/2500=4(47/100)².             (11)

No such radius is a strict full witness. Radii r>=1/2 are excluded by the original cap(1): apply it to the closed square of side2r containing the open receiving square. Therefore

    E intersect H(sigma+nu) is empty,
    E subset M.                                    (12)

The proxy has replaced a potentially tiny captured old prefix by all100 local old atoms. Equation(12) proves a genuine positive-area overestimate. It neither calculates the complete actual deletion nor disproves any bound on that deletion.

## 5. The collapse is not a feasible replacement of the full grid

There are two different objects: the original full input satisfying(1), and the proxy used to upper-envelope some of its mixed witnesses. They must not be identified. If all100 old points in P are actually replaced by mass1/4 at the origin while retaining the old grid outside P, the resulting synthetic full input violates a unit-square cap.

Take

    Q=[-71/200,129/200] x [-171/200,29/200].          (13)

The original Q contains13 columns and3 rows of old atoms, plus both outside atoms, so its mass is39/400+3/4=339/400. Among those39 old atoms,30 belong to P and9 are outside P. After the replacement, Q contains the collapsed mass1/4, the9 retained old atoms and outside mass3/4, giving

    1/4+9/400+3/4=409/400>1.                        (14)

Thus this result retains full caps as hypotheses on the original input; it does not assert that the synthetic full replacement preserves them. The isolated collapsed three-atom measure has total mass1 and hence satisfies caps with normalization1, but its selected old mass is1/4 and is not a maximum eligible mass. That separate observation does not repair the proposed half-local-mass payment.

## 6. Executed checks and remaining mathematical task

The production component computes the complete three mixed receivers and all three residual receivers by two directed72-bit union algorithms, obtaining an interval consistent with(3). It also checks the original n20 candidate a=31/100,d=83/100,b1=1/8,b2=1/2, whose proxy area is approximately0.110722446721674<1/8. That failed candidate remains a failed test; the corrected rational input above is separately identified. The uncapped mixed-only obstruction belongs to the companion proof and its own calculation.

The present author's `computation/round077/check_qualified_proxy_scalars.py` was actually run once using Python3.15.0 Fraction. It checks the squared orderings, expansion coefficients, rational strict lower bound, receiver-cell margins, all400 old positions for the unit-window count, and the connecting-edge/nonterminal inequalities. PID83088 exited0 and was confirmed absent. It does not repeat the440 cap rows,45302 closure enumeration,100 activity rows or two complete area algorithms. The input, result and receipt are saved beside the script.

The actual uniform-grid corner payment, the joint-active two-outside case and the exact centered square weak(1,1) constant remain open. The practical restriction established here is precise: even after requiring a legal original full grid and keeping Gamma, r<1/2 and the true residual, one cannot replace every captured internal prefix by the total m_P and then pay its unified mixed proxy by m_P/2. A successful active-branch proof must retain more original prefix geometry or a different payment mechanism.
