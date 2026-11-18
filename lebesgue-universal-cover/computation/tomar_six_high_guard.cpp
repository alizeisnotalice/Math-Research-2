#include "fan-exposed.hpp"

struct DynamicGuard {
    Tester& t;
    const L* directions;
    const L* disk;
    const L* endpoints;
    const L* pointers;
    const L* triples;
    Box root_lo,root_hi;
    int n, cells, depth_limit, calls_limit;
    L calls=0, old=0, added=0, depth=0, tried=0;
    Box unknown_lo{},unknown_hi{};
    int reason=0;
    std::vector<L> needed;

    static constexpr L U=L(1)<<22;
    static constexpr L Scale=2*Q*Q;
    static L up_signed(W p,L q){
        return (p>=0 ? (p+q-1)/q : p/q);
    }

    bool excluded(const Box& lo,const Box& hi){
        W e=t.goal-t.last;
        if(e<=0)return true;
        if(e>W(Scale)/4000)return false;
        ++tried;
        needed.assign(n,-1);
        auto threshold=[&](int j){
            if(needed[j]>=0)return needed[j];
            L b=disk[j];
            L x=directions[2*j],y=directions[2*j+1];
            for(int i=0;i<4;++i){
                L h=endpoints[(2*i)*n+j],p=endpoints[(2*i+1)*n+j];
                if(i){
                    int k=3*i-1;
                    L q=p>=h ? hi[k] : lo[k];
                    L d=root_hi[k]-root_lo[k];
                    if(q<root_lo[k] || q>root_hi[k] || d<=0)std::abort();
                    h=up_signed(W(h)*(root_hi[k]-q)+W(p)*(q-root_lo[k]),d);
                }
                int k=3*i;
                W pmax=W(x)*(x>=0?hi[k]:lo[k])+W(y)*(y>=0?hi[k+1]:lo[k+1]);
                h+=up_signed(pmax, (Fine/Q)*U);
                b=std::max(b,h);
            }
            // Area is tested at a common dyadic cap scale C.  The upward
            // rounded b and eta make the comparison stronger than needed.
            constexpr L C=L(1)<<20, C2=C*C;
            L B=up_signed(W(b)*U,(U-2)*(Q/C));
            L area=up_signed(e*C2,Scale);
            auto sufficient=[&](L mid){
                L sum=B+mid;
                W lhs=W(C2)*mid*mid*mid*mid;
                W rhs=W(area)*area*(W(4)*sum*sum-C2);
                return lhs>=rhs;
            };
            // Floating evaluation merely supplies a guess. Repair it by exact tests.
            double bound=double(B)/C, gap=double(area)/double(C2), guess=.020;
            for(int i=0;i<7;++i)guess=std::sqrt(2*gap*std::sqrt((bound+guess)*(bound+guess)-.25));
            L r=std::max<L>(1,L(std::ceil(guess*C)));
            while(!sufficient(r))++r;
            while(r>1 && sufficient(r-1))--r;
            L d=r*(Q/C);
            return needed[j]=b+up_signed(W(d)*(U+2),U);
        };
        for(int i=0;i<cells;++i){
            bool ok=false;
            for(L k=pointers[i];k<pointers[i+1];++k){
                const L* q=triples+8*k;
                W support=W(q[3])*threshold(q[0])+W(q[4])*threshold(q[1])+W(q[5])*threshold(q[2]);
                if(W(q[6])*q[7]>support){ok=true;break;}
            }
            if(!ok)return false;
        }
        return true;
    }

    bool stop(const Box& lo,const Box& hi,int why){
        unknown_lo=lo;unknown_hi=hi;reason=why;
        return false;
    }

    bool walk(Box& lo,Box& hi,int level){
        if(calls==calls_limit)return stop(lo,hi,1);
        ++calls; depth=std::max(depth,L(level));
        int ret=t.check(lo,hi);
        if(ret<0)return stop(lo,hi,4);
        if(ret==1){++old;return true;}
        if(excluded(lo,hi)){++added;return true;}
        if(level==depth_limit)return stop(lo,hi,2);
        int k=t.split_axis;
        L a=lo[k],b=hi[k];
        if(b<=a || (a+b)%2)return stop(lo,hi,3);
        L mid=(a+b)/2;
        hi[k]=mid;bool left=walk(lo,hi,level+1);hi[k]=b;
        if(!left)return false;
        lo[k]=mid;bool right=walk(lo,hi,level+1);lo[k]=a;
        return right;
    }
};

extern "C" int dynamic_guard_run(const L* lower,const L* upper,const L* original,
  const std::int32_t* table,const std::int32_t* counts,L goal,int depth,int cap,
  const L* dirs,const L* disk,const L* endpoints,const L* pointers,const L* triples,
  int normals,int cells,L* out){
    Tester tester(table,counts,original,goal);
    Box lo,hi;std::copy(lower,lower+11,lo.begin());std::copy(upper,upper+11,hi.begin());
    DynamicGuard run{tester,dirs,disk,endpoints,pointers,triples,lo,hi,normals,cells,depth,cap};
    bool closed=run.walk(lo,hi,0);
    out[0]=closed;out[1]=run.calls;out[2]=run.old;out[3]=run.added;out[4]=run.depth;out[5]=run.reason;out[6]=run.tried;
    if(!closed){std::copy(run.unknown_lo.begin(),run.unknown_lo.end(),out+7);std::copy(run.unknown_hi.begin(),run.unknown_hi.end(),out+18);}
    return 0;
}
