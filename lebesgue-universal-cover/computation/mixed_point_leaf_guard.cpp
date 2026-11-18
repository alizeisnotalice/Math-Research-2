#include "tomar_six_083730_guard.cpp"

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

    bool fail(const Box& lo,const Box& hi,int code) {
        failure=code;failure_lo=lo;failure_hi=hi;return false;
    }
    void emit(const Box& lo,const Box& hi) {
        L* dest=trace+22*scanned++;
        std::copy(lo.begin(),lo.end(),dest);
        std::copy(hi.begin(),hi.end(),dest+11);
    }

    bool walk(Box& lo,Box& hi,int level){
        if(nodes==budget)return fail(lo,hi,1);
        ++nodes;
        depth=std::max(depth,L(level));
        int found=old.check(lo,hi);
        if(found<0)return fail(lo,hi,4);
        if(found==1){++old_stops;return true;}
        if(p.excluded(lo,hi)){++p_stops;return true;}
        if(mode&2) {
            const int g=bank.covers(lo,hi);
            if(g>=0){++mixed_stops;++uses[g];++geometric[bank.data[13*g+12]];return true;}
        }
        if((mode&1) && old.goal-old.last<=W(DynamicGuard::Scale)/25000 && q.excluded(lo,hi)){
            ++q_stops;
            if(trace_kind==0 && scanned<trace_capacity && q_stops%stride==0)emit(lo,hi);
            return true;
        }
        if(trace_kind==1 && scanned<trace_capacity && level>=45) {
            std::array<L,12> mid;
            for(int j=0;j<11;++j)mid[j]=lo[j]+hi[j];
            mid[11]=(L(1)<<45)/5000;
            if(CenterBank::contains(mid.data(),lo,hi) && ++internal_samples%stride==0)emit(lo,hi);
        }
        if(level==max_depth)return fail(lo,hi,2);
        const int k=old.split_axis;
        const L a=lo[k],b=hi[k];
        if(b<=a || (a+b)%2)return fail(lo,hi,3);
        const L m=(a+b)/2;
        hi[k]=m;
        const bool left=walk(lo,hi,level+1);
        hi[k]=b;
        if(!left)return false;
        lo[k]=m;
        const bool right=walk(lo,hi,level+1);
        lo[k]=a;
        return right;
    }
};

extern "C" int mixed_point_leaf_run(const L* lower,const L* upper,const L* original,
  const std::int32_t* table,const std::int32_t* counts,L goal,int depth,int cap,int mode,
  const L* p_dirs,const L* p_disk,const L* p_endpoints,const L* p_pointers,const L* p_triples,int pn,
  const L* q_dirs,const L* q_disk,const L* q_endpoints,const L* q_pointers,const L* q_triples,int qn,
  const L* centers,int center_count,L* trace,int trace_capacity,int stride,int trace_kind,L* out){
    Tester tester(table,counts,original,goal);
    Box lo,hi;
    std::copy(lower,lower+11,lo.begin());
    std::copy(upper,upper+11,hi.begin());
    DynamicGuard p{tester,p_dirs,p_disk,p_endpoints,p_pointers,p_triples,lo,hi,pn,384,depth,cap,2};
    DynamicGuard q{tester,q_dirs,q_disk,q_endpoints,q_pointers,q_triples,lo,hi,qn,384,depth,cap,2};
    p.init_collective();
    q.init_collective();
    MixedSearch search{tester,p,q,{centers,center_count},cap,depth,mode,trace_capacity,stride,trace_kind,trace,out+35};
    const bool pass=search.walk(lo,hi,0);
    L tally[]={pass,search.nodes,search.old_stops,search.p_stops,search.q_stops,search.mixed_stops,
              search.depth,search.failure,p.tried,q.tried,search.scanned};
    std::copy(tally,tally+11,out);
    std::copy(search.geometric,search.geometric+2,out+33);
    if(!pass) {
        std::copy(search.failure_lo.begin(),search.failure_lo.end(),out+11);
        std::copy(search.failure_hi.begin(),search.failure_hi.end(),out+22);
    }
    return 0;
}
