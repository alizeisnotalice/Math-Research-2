#include "tomar_six_083730_guard.cpp"

struct PairSearch {
    Tester& old;
    DynamicGuard& p;
    DynamicGuard& q;
    int budget, max_depth;
    bool use_q;
    L nodes=0, old_stops=0, p_stops=0, q_stops=0, depth=0;
    Box failure_lo{}, failure_hi{};
    int failure_code=0;

    bool fail(const Box& lo,const Box& hi,int why){
        failure_lo=lo;
        failure_hi=hi;
        failure_code=why;
        return false;
    }

    bool walk(Box& lo,Box& hi,int level){
        if(nodes==budget) return fail(lo,hi,1);
        ++nodes;
        depth=std::max(depth,L(level));
        int status=old.check(lo,hi);
        if(status<0) return fail(lo,hi,4);
        if(status==1){++old_stops;return true;}
        if(p.excluded(lo,hi)){++p_stops;return true;}
        if(use_q && old.goal-old.last<=W(DynamicGuard::Scale)/25000 && q.excluded(lo,hi)){
            ++q_stops;
            return true;
        }
        if(level==max_depth) return fail(lo,hi,2);
        int k=old.split_axis;
        L a=lo[k],b=hi[k];
        if(b<=a || (a+b)%2) return fail(lo,hi,3);
        L mid=(a+b)/2;
        hi[k]=mid;
        bool left=walk(lo,hi,level+1);
        hi[k]=b;
        if(!left) return false;
        lo[k]=mid;
        bool right=walk(lo,hi,level+1);
        lo[k]=a;
        return right;
    }
};

extern "C" int dual_point_guard_run(const L* lower,const L* upper,const L* original,
  const std::int32_t* table,const std::int32_t* counts,L goal,int depth,int cap,int use_q,
  const L* p_dirs,const L* p_disk,const L* p_endpoints,const L* p_pointers,const L* p_triples,int pn,
  const L* q_dirs,const L* q_disk,const L* q_endpoints,const L* q_pointers,const L* q_triples,int qn,
  L* out){
    Tester tester(table,counts,original,goal);
    Box lo,hi;
    std::copy(lower,lower+11,lo.begin());
    std::copy(upper,upper+11,hi.begin());
    DynamicGuard p{tester,p_dirs,p_disk,p_endpoints,p_pointers,p_triples,lo,hi,pn,384,depth,cap,2};
    DynamicGuard q{tester,q_dirs,q_disk,q_endpoints,q_pointers,q_triples,lo,hi,qn,384,depth,cap,2};
    p.init_collective();
    q.init_collective();
    PairSearch run{tester,p,q,cap,depth,use_q>0};
    bool pass=run.walk(lo,hi,0);
    out[0]=pass;out[1]=run.nodes;out[2]=run.old_stops;out[3]=run.p_stops;out[4]=run.q_stops;
    out[5]=run.depth;out[6]=run.failure_code;out[7]=p.tried;out[8]=q.tried;
    out[9]=p.collective_boxes;out[10]=q.collective_boxes;
    if(!pass){
        std::copy(run.failure_lo.begin(),run.failure_lo.end(),out+11);
        std::copy(run.failure_hi.begin(),run.failure_hi.end(),out+22);
    }
    return 0;
}
