#include <cstdio>
#include <functional>
#include "angle_mesh_dynamic.hpp"

using Reanchor = int (*)(const L*, const L*, L*, L*);
using CenterProof = L (*)(const L*, const L*);

struct GuardAnchor {
    DynamicGuard& guard;
    const L* endpoints=guard.endpoints;
    Box lo=guard.root_lo, hi=guard.root_hi;
    AngleMesh mesh=guard.mesh;

    void reset(const L* fresh,const Box& lower,const Box& upper) {
        guard.endpoints=fresh;
        guard.root_lo=lower;
        guard.root_hi=upper;
        guard.mesh={};
    }
    ~GuardAnchor() {
        guard.endpoints=endpoints;
        guard.root_lo=lo;
        guard.root_hi=hi;
        guard.mesh=mesh;
    }
};

struct CenterBank {
    const L* data;
    int count;
    static bool contains(const L* row,const Box& lo,const Box& hi) {
        const L radius=row[11];
        const auto dist=[&](int j) {
            return std::max(std::llabs(2*lo[j]-row[j]), std::llabs(2*hi[j]-row[j]));
        };
        constexpr int slope[4]={0,163,115,89};
        for(int b=0;b<4;++b) {
            const int i=3*b;
            const L x=dist(i),y=dist(i+1);
            if(std::max(x,y)>radius)return false;
            const L rotation=b ? slope[b]*dist(i-1) : 0;
            const L allowance=1000*radius-rotation;
            if(allowance<0)return false;
            if(W(1000)*1000*(W(x)*x+W(y)*y)>W(allowance)*allowance)return false;
        }
        return true;
    }
    int covers(const Box& lo,const Box& hi) const {
        for(int j=0;j<count;++j)
            if(contains(data+13*j,lo,hi))return j;
        return -1;
    }
};

struct MixedSearch {
    Tester& old;
    DynamicGuard& p;
    DynamicGuard& q;
    CenterBank bank;
    int budget, max_depth, mode, trace_capacity, stride, trace_kind;
    L* trace;
    L* uses;
    L internal_samples=0;
    L nodes=0, old_stops=0, p_stops=0, q_stops=0, mixed_stops=0, depth=0, scanned=0;
    L geometric[2]{};
    int failure=0;
    Box failure_lo{},failure_hi{};
    L* frontier;
    int reanchor_level, reanchor_stride;
    Reanchor reanchor;
    CenterProof center;
    int center_level, center_limit, sweep_limit;
    std::vector<L> path;
    std::vector<std::array<L,13>> new_centers;
    L pending=0, anchors=0;
    L pending_visited=0;
    L center_calls=0, center_stops=0;

    void keep(const Box& lo,const Box& hi,int level) {
        L* row=frontier + pending++*(23+max_depth);
        row[0]=level;
        std::copy(lo.begin(),lo.end(),row+1);
        std::copy(hi.begin(),hi.end(),row+12);
        std::copy(path.begin(),path.begin()+level,row+23);
    }
    bool fail(const Box& lo,const Box& hi,int level,int code) {
        keep(lo,hi,level);
        if(code!=1)++pending_visited;
        if(!failure){
            failure=code;failure_lo=lo;failure_hi=hi;
        }
        return (code==2 || code==3) && pending<sweep_limit;
    }
    void emit(const Box& lo,const Box& hi) {
        L* dest=trace+22*scanned++;
        std::copy(lo.begin(),lo.end(),dest);
        std::copy(hi.begin(),hi.end(),dest+11);
    }

    bool guard(const Box& lo,const Box& hi) {
        if(!(mode&4) && p.excluded(lo,hi)){++p_stops;return true;}
        if(mode&2) {
            const int g=bank.covers(lo,hi);
            if(g>=0){++mixed_stops;++uses[g];++geometric[bank.data[13*g+12]];return true;}
        }
        for(const auto& row:new_centers) {
            if(CenterBank::contains(row.data(),lo,hi)) {
                ++mixed_stops;++center_stops;++geometric[0];
                return true;
            }
        }
        if((mode&1) && old.goal-old.last<=W(DynamicGuard::Scale)/25000 && q.excluded(lo,hi)){
            ++q_stops;
            if(trace_kind==0 && scanned<trace_capacity && q_stops%stride==0)emit(lo,hi);
            return true;
        }
        return false;
    }

    bool walk(Box& lo,Box& hi,int level,int next_reanchor) {
        if(nodes==budget)return fail(lo,hi,level,1);
        ++nodes;
        if(nodes%500000==0) {
            std::fprintf(stderr,
                "mode=%d nodes=%lld five=%lld P30=%lld Q11=%lld depth=%lld anchors=%lld\n",
                mode,(long long)nodes,(long long)old_stops,(long long)p_stops,
                (long long)q_stops,(long long)depth,(long long)anchors);
        }
        depth=std::max(depth,L(level));
        const int found=old.check(lo,hi);
        if(found<0)return fail(lo,hi,level,4);
        if(found==1){++old_stops;return true;}
        if(guard(lo,hi))return true;
        if(reanchor_level>0 && !(mode&4) && level>=next_reanchor && level<max_depth) {
            std::vector<L> p_end(8*p.n),q_end((mode&1)?8*q.n:0);
            if(reanchor(lo.data(),hi.data(),p_end.data(),q_end.data()))
                return fail(lo,hi,level,5);
            ++anchors;
            GuardAnchor p_scope{p},q_scope{q};
            p_scope.reset(p_end.data(),lo,hi);
            if(mode&1)q_scope.reset(q_end.data(),lo,hi);
            if(guard(lo,hi))return true;
            return split(lo,hi,level,level+reanchor_stride);
        }
        return split(lo,hi,level,next_reanchor);
    }

    bool split(Box& lo,Box& hi,int level,int next_reanchor) {
        if(trace_kind==1 && scanned<trace_capacity && level>=45) {
            std::array<L,12> mid;
            for(int j=0;j<11;++j)mid[j]=lo[j]+hi[j];
            mid[11]=(L(1)<<45)/5000;
            if(CenterBank::contains(mid.data(),lo,hi) && ++internal_samples%stride==0)emit(lo,hi);
        }
        const int k=old.split_axis;
        const L a=lo[k],b=hi[k];
        const bool finest=b<=a || (a+b)%2;
        if(center_limit>center_calls && !(mode&4) && (level>=center_level || finest)) {
            ++center_calls;
            const L radius=center(lo.data(),hi.data());
            if(radius<0)return fail(lo,hi,level,6);
            if(radius>0) {
                std::array<L,13> row{};
                std::transform(lo.begin(),lo.end(),hi.begin(),row.begin(),std::plus<L>{});
                row[11]=radius;
                new_centers.push_back(row);
                if(CenterBank::contains(row.data(),lo,hi)) {
                    ++mixed_stops;++center_stops;++geometric[0];
                    return true;
                }
            }
        }
        if(level==max_depth)return fail(lo,hi,level,2);
        if(finest)return fail(lo,hi,level,3);
        const L m=(a+b)/2;
        hi[k]=m;
        path[level]=2*k;
        const bool left=walk(lo,hi,level+1,next_reanchor);
        hi[k]=b;
        lo[k]=m;
        path[level]=2*k+1;
        if(!left) {
            keep(lo,hi,level+1);
            lo[k]=a;
            return false;
        }
        const bool right=walk(lo,hi,level+1,next_reanchor);
        lo[k]=a;
        return right;
    }
};

extern "C" int mixed_point_leaf_run(const L* lower,const L* upper,const L* original,
  const std::int32_t* table,const std::int32_t* counts,L goal,int depth,int cap,int mode,
  const L* p_dirs,const L* p_disk,const L* p_endpoints,const L* p_pointers,const L* p_triples,int pn,
  const L* q_dirs,const L* q_disk,const L* q_endpoints,const L* q_pointers,const L* q_triples,int qn,
  const L* centers,int center_count,L* trace,int trace_capacity,int stride,int trace_kind,L* out,
  const L* knots,const L* p_values,const L* q_values,int panels,int q_mesh,
  L* frontier,Reanchor reanchor,int reanchor_level,int reanchor_stride,
  CenterProof center,int center_level,int center_limit,int sweep_limit){
    Tester tester(table,counts,original,goal);
    Box lo,hi;
    std::copy(lower,lower+11,lo.begin());
    std::copy(upper,upper+11,hi.begin());
    DynamicGuard p{tester,p_dirs,p_disk,p_endpoints,p_pointers,p_triples,lo,hi,pn,384,depth,cap,2};
    DynamicGuard q{tester,q_dirs,q_disk,q_endpoints,q_pointers,q_triples,lo,hi,qn,384,depth,cap,2};
    if(panels) {
        p.mesh={knots,p_values,panels};
        if(q_mesh && (mode&1))q.mesh={knots,q_values,panels};
    }
    p.init_collective();
    q.init_collective();
    MixedSearch search{tester,p,q,{centers,center_count},cap,depth,mode,trace_capacity,stride,trace_kind,trace,out+35};
    search.frontier=frontier;
    search.path.resize(depth);
    search.reanchor=reanchor;
    search.reanchor_level=reanchor_level;
    search.reanchor_stride=reanchor_stride;
    search.center=center;
    search.center_level=center_level;
    search.center_limit=center_limit;
    search.sweep_limit=sweep_limit;
    const bool pass=search.walk(lo,hi,0,reanchor_level) && search.pending==0;
    L tally[]={pass,search.nodes,search.old_stops,search.p_stops,search.q_stops,search.mixed_stops,
              search.depth,search.failure,p.tried,q.tried,search.scanned};
    std::copy(tally,tally+11,out);
    std::copy(search.geometric,search.geometric+2,out+33);
    out[35+center_count]=search.pending;
    out[36+center_count]=search.anchors;
    out[37+center_count]=search.center_calls;
    out[38+center_count]=search.center_stops;
    out[39+center_count]=p.edge_boxes;
    out[40+center_count]=q.edge_boxes;
    out[41+center_count]=search.pending_visited;
    if(!pass) {
        std::copy(search.failure_lo.begin(),search.failure_lo.end(),out+11);
        std::copy(search.failure_hi.begin(),search.failure_hi.end(),out+22);
    }
    return 0;
}
