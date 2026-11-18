// Whole-box chord caps in the actual five-body hull, using the already audited
// moved witnesses of Tester. Heights and tangential coordinates are in Q units
// along vectors scaled by U, before normalization.
struct LongEdgeSupport {
    const Tester& test;
    static constexpr L U=L(1)<<22;
    static_assert(W(U+2)*(U+2)*1023 < W(U)*U*1024);
    constexpr static std::array<L,6> band{Q/2500,Q/1000,Q/500,Q/250,Q/125,2*Q/125};
    std::vector<L> displacement;

    static L upper(W n,L d) {
        return L(n>=0 ? (n+d-1)/d : n/d);
    }

    L bound(L vx,L vy,L high,W gap,L from_disk) {
        std::array<L,6> left,right;
        left.fill(std::numeric_limits<L>::max()/4);
        right.fill(std::numeric_limits<L>::min()/4);
        if(displacement.empty()) {
            for(const auto& p:test.hull)
                displacement.push_back(upper(W(test.perturb(p))*(U+2),U));
        }
        for(std::size_t k=0;k<test.hull.size();++k) {
            const auto& p=test.hull[k];
            const L error=displacement[k];
            const W n=W(vx)*p.x+W(vy)*p.y;
            const L hlo=-upper(-n,U)-error;
            const L hhi=upper(n,U)+error;
            if(hlo<high-band.back()) continue;
            const W v=W(vx)*p.y-W(vy)*p.x;
            const L tlo=-upper(-v,U)-error, thi=upper(v,U)+error;
            for(std::size_t j=0;j<band.size();++j) {
                const L h=high-band[j];
                if(h<Q*45/100 || hlo<h) continue;
                // Scale each point radially to height h. Its change in
                // tangential projection is bounded by 2*(hhi-h) since the
                // radius-0.707 bound gives |t| <= 2*h in any admissible host.
                const L loss=2*(hhi-h);
                left[j]=std::min(left[j],thi+loss);
                right[j]=std::max(right[j],tlo-loss);
            }
        }
        L best=from_disk;
        for(std::size_t j=0;j<band.size();++j) {
            const L length=right[j]-left[j];
            if(length<=0) continue;
            const W c=W(1023)*length, a=W(1024)*gap;
            const auto enough=[&](L r){return c*r*r>=a*(W(r)+band[j]);};
            if(best<=1 || !enough(best-1)) continue;
            const double A=double(a), C=double(c);
            L r=std::max<L>(1,L(std::ceil((A+std::sqrt(A*A+4*A*C*band[j]))/(2*C))));
            while(!enough(r))++r;
            while(r>1 && enough(r-1))--r;
            best=std::min(best,r);
        }
        return best;
    }
};
