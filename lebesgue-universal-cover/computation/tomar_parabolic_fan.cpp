// Integer fan-area lower bound on whole 11D dyadic boxes; coordinates in units 2^-28.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <vector>
using L=std::int64_t;using W=__int128_t;
constexpr L Q=L(1)<<28,Fine=L(1)<<44,ratio=Fine/Q;
struct P { L x,y,a,b; int body; };
struct V { L x,y; };
using Box=std::array<L,11>;
W det(L ax,L ay,L bx,L by){return W(ax)*by-W(ay)*bx;}
W area(const P& a,const P& b,const P& c){return det(b.x-a.x,b.y-a.y,c.x-a.x,c.y-a.y);}
L norm(L x,L y){return std::abs(x)+std::abs(y);}
W up(W a,L b){
    if(a<0 || b<=0) std::abort();
    return (a+b-1)/b;
}
L small(W a){
    if(a<0 || a>std::numeric_limits<L>::max()) std::abort();
    return static_cast<L>(a);
}
W ab(W x){return x<0?-x:x;}
struct Tester {
    const std::int32_t* table;const std::int32_t* sizes;const L* original_hi;
    L goal;int total=0;int split_axis=0;
    std::vector<P> pts,hull,conv;
    L hx[5]{},hy[5]{},theta[5]{},tn[5]{},tw[5]{},rad[5]{};
    Tester(const std::int32_t* p,const std::int32_t* n,const L* h,L t):table(p),sizes(n),original_hi(h),goal(t){
        for(int k=0;k<5;++k){if(sizes[k]<3 || sizes[k]>1024)std::abort();total+=sizes[k];}
        if(total>4096)std::abort();
        for(int k=0;k<11;++k)if(h[k]<=-Fine || h[k]>=Fine)std::abort();
        for(int k: {2,5,8})if(h[k]<=0)std::abort();
        pts.reserve(total);hull.reserve(2*total);conv.reserve(total);
    }
    L perturb(const P& p) const {int b=p.body;return small(W(tw[b])+up(W(norm(p.a,p.b))*theta[b],Q));}
    bool oriented(const P& a,const P& b) const {
        W z=det(a.x,a.y,b.x,b.y);
        if(z<=0)return false;
        if(a.body==b.body){
            int i=a.body;L chord=norm(a.a-b.a,a.b-b.b);
            return z>W(tw[i])*chord+W(tn[i])*up(W(chord)*theta[i],Q);
        }
        L da=perturb(a),db=perturb(b);
        return z>W(norm(a.x,a.y))*db+W(norm(b.x,b.y))*da+W(da)*db;
    }
    W curved_gain(const P& a,const P& b) const {
        if(a.body!=b.body)return 0;
        W cross=det(a.x,a.y,b.x,b.y);
        if(!a.body){
            // The actual radius-one-half disk contains the entire radial sector.
            // For 0<alpha<pi, alpha >= sin(alpha)+sin(alpha)^3/6.
            // Both fractional radii are below 1/2, so 4*cross/Q^2 <= sin(alpha).
            W sine=(4*cross)/Q;
            return (sine*sine*sine)/(24*Q);
        }
        L vx=b.a-a.a,vy=b.b-a.b,chord=norm(vx,vy);
        if(!chord)return 0;
        W length_sq=W(vx)*vx+W(vy)*vy;
        // The floating guess is used only to recover the exact integer square root.
        L length=static_cast<L>(std::sqrt(static_cast<double>(length_sq)));
        while(W(length)*length>length_sq)--length;
        while(W(length+1)*(length+1)<=length_sq)++length;
        int body=a.body;
        W zlo=cross-W(tw[body])*chord-W(tn[body])*up(W(chord)*theta[body],Q);
        if(zlo<=0)return 0;
        L radius_limit=std::max(norm(a.x,a.y)+perturb(a),norm(b.x,b.y)+perturb(b));
        if(!radius_limit)return 0;
        // The strip x_s-u*J(b-a), 0<=u<=kappa*s*(1-s), lies in all defining unit disks.
        // The determinant bound places the whole parabolic strip in this angular wedge.
        W kappa=std::min<W>(length/2,(zlo*Q)/(W(radius_limit)*chord));
        return (kappa*length_sq)/(3*Q);
    }
    int check(const Box& lo,const Box& hi){
        pts.clear();hull.clear();conv.clear();int off=0;
        for(int b=0;b<5;++b){
            int idx=b?3*(b-1):0;L tx=0,ty=0;int row=0;
            if(b){
                tx=(lo[idx]+hi[idx])/(2*ratio);ty=(lo[idx+1]+hi[idx+1])/(2*ratio);
                // Nominal translations are rounded to Q. Enlarge their half-widths exactly.
                hx[b]=small(up(std::max(std::abs(lo[idx]-tx*ratio),std::abs(hi[idx]-tx*ratio)),ratio));
                hy[b]=small(up(std::max(std::abs(lo[idx+1]-ty*ratio),std::abs(hi[idx+1]-ty*ratio)),ratio));
                if(b>1){
                    int k=idx-1;L num=(lo[k]+hi[k])*512;
                    row=int((num+original_hi[k]/2)/original_hi[k]);if(row<0 || row>1024)return -1;
                    L anchor=L(row)*original_hi[k];
                    theta[b]=small(up(std::max(std::abs(1024*lo[k]-anchor),std::abs(1024*hi[k]-anchor)),1024*ratio));
                }
            }
            tw[b]=hx[b]+hy[b];tn[b]=norm(tx,ty);
            const auto* block=table+2*(row*total+off);
            for(int j=0;j<sizes[b];++j){
                L a=block[2*j],c=block[2*j+1];
                pts.push_back(P{tx+a,ty+c,a,c,b});
            }
            off+=sizes[b];
        }
        std::sort(pts.begin(),pts.end(),[](const P&a,const P&b){return a.x!=b.x?a.x<b.x:(a.y!=b.y?a.y<b.y:a.body<b.body);});
        pts.erase(std::unique(pts.begin(),pts.end(),[](const P&a,const P&b){return a.x==b.x && a.y==b.y;}),pts.end());
        for(auto& p:pts){
            while(hull.size()>=2 && area(hull[hull.size()-2],hull.back(),p)<=0)hull.pop_back();
            hull.push_back(p);
        }
        auto lower=hull.size();
        for(int i=int(pts.size())-2;i>=0;--i){
            while(hull.size()>lower && area(hull[hull.size()-2],hull.back(),pts[i])<=0)hull.pop_back();
            hull.push_back(pts[i]);
        }
        hull.pop_back();conv=hull;
        // Delete a candidate only when it is cheaper than the other endpoint of the bad edge.
        bool done=false;
        while(!done && conv.size()>=3){
            done=true;auto n=conv.size();
            for(unsigned i=0;i<n;++i){
                unsigned j=(i+1)%n;
                if(oriented(conv[i],conv[j]))continue;
                auto pre=(i+n-1)%n,post=(j+1)%n;
                auto loss_i=area(conv[pre],conv[i],conv[j]),loss_j=area(conv[i],conv[j],conv[post]);
                conv.erase(conv.begin()+(loss_i<loss_j?i:j));done=false;break;
            }
        }
        if(conv.size()<3){
            split_axis=0;for(int k=1;k<11;++k)if(hi[k]-lo[k]>hi[split_axis]-lo[split_axis])split_axis=k;
            return 0;
        }
        W twice=0,error=0,grot[5]{};V g[5]{};
        auto n=conv.size();
        for(unsigned i=0;i<n;++i){
            const P& a=conv[i];const P& next=conv[(i+1)%n];const P& prev=conv[(i+n-1)%n];
            L dx=next.x-prev.x,dy=next.y-prev.y;
            twice+=det(a.x,a.y,next.x,next.y);
            twice+=curved_gain(a,next);
            int k=a.body;g[k].x+=dx;g[k].y+=dy;
            grot[k]-=W(a.a)*dx+W(a.b)*dy;
            L d1=small(up(W(norm(a.a,a.b))*theta[k],Q));
            error+=W(up(W(d1)*theta[k],2*Q))*norm(dx,dy);
            if(k!=next.body){
                error+=W(perturb(a))*perturb(next);
            }else{
                L chord=norm(next.a-a.a,next.b-a.b);
                error+=W(tw[k])*up(W(chord)*theta[k],Q);
                W signed_cross=det(a.a,a.b,next.a,next.b);
                if(signed_cross<0){
                    auto first=up((-signed_cross)*theta[k],Q);
                    error+=up(W(first)*theta[k],Q);
                }
            }
        }
        for(int k=0;k<5;++k){
            error+=W(std::abs(g[k].y))*hx[k]+W(std::abs(g[k].x))*hy[k];
            error+=up(ab(grot[k])*theta[k],Q);
        }
        W best=-1;
        for(int b=1;b<5;++b){
            int k=3*(b-1);
            auto consider=[&](int index,W score){if(score>best){split_axis=index;best=score;}};
            consider(k,W(hx[b])*(std::abs(g[b].y)+Q/30));
            consider(k+1,W(hy[b])*(std::abs(g[b].x)+Q/30));
            if(b>1)consider(k-1,W(theta[b])*(ab(grot[b])/Q+Q/100));
        }
        return twice-error>=goal ? 1 : 0;
    }
};
struct Runner {
    Tester& tester;const std::uint8_t* codes;const std::vector<int>& stop;
    const std::vector<int>& active;std::uint8_t* mask;L* stats;
    bool walk(int i,Box& lo,Box& hi,int depth){
        if(depth>256)return false;
        int count=active[stop[i]]-active[i];if(!count)return true;
        int k=codes[i];
        if(k>=253 || (count>=32 && depth%5==0)){
            ++stats[3];int ret=tester.check(lo,hi);if(ret<0)return false;
            if(ret>0){stats[0]+=count;++stats[1];return true;}
        }
        if(k>=253){mask[i]=1;++stats[2];return true;}
        L l=lo[k],h=hi[k];if((l+h)%2)return false;
        L mid=(l+h)/2;hi[k]=mid;
        if(!walk(i+1,lo,hi,depth+1))return false;
        hi[k]=h;lo[k]=mid;
        if(!walk(stop[i+1],lo,hi,depth+1))return false;
        lo[k]=l;return true;
    }
};
extern "C" int fan_trim(const std::uint8_t* codes,int n,const L* lower,const L* upper,
    const L* original_hi,const std::uint8_t* old,std::uint8_t* next,L* stats,
    const std::int32_t* table,const std::int32_t* counts,L goal){
    if(n<1)return 1;
    std::vector<int> stop(n),active(n+1);
    for(int i=0;i<n;++i){if(old[i]&&codes[i]<253)return 2;active[i+1]=active[i]+(old[i]?1:0);}
    for(int i=n-1;i>=0;--i){
        if(codes[i]>=253)stop[i]=i+1;
        else{if(codes[i]>=11 || i+1>=n || stop[i+1]>=n)return 3;stop[i]=stop[stop[i+1]];}
    }
    if(stop[0]!=n)return 4;
    Tester tester(table,counts,original_hi,goal);
    Runner runner{tester,codes,stop,active,next,stats};
    Box lo,hi;std::copy(lower,lower+11,lo.begin());std::copy(upper,upper+11,hi.begin());
    return runner.walk(0,lo,hi,0)?0:5;
}

bool split_probe(Tester& t,Box& lo,Box& hi,int depth,int maxd,int maxcalls,L* stat){
    if(stat[0]>=maxcalls)return false;
    ++stat[0];stat[1]=std::max(stat[1],L(depth));
    if(t.check(lo,hi)==1){++stat[2];return true;}
    if(depth==maxd)return false;
    int k=t.split_axis;auto l=lo[k],h=hi[k];if((l+h)%2 || h<=l)return false;
    auto mid=(l+h)/2;hi[k]=mid;
    bool ok=split_probe(t,lo,hi,depth+1,maxd,maxcalls,stat);hi[k]=h;
    if(!ok)return false;
    lo[k]=mid;ok=split_probe(t,lo,hi,depth+1,maxd,maxcalls,stat);lo[k]=l;
    return ok;
}
extern "C" int fan_probe(int n,const L* lo,const L* hi,const L* original_hi,const std::int32_t* table,
    const std::int32_t* counts,L goal,int depth,int calls,L* results){
    Tester t(table,counts,original_hi,goal);
    for(int i=0;i<n;++i){
        Box l,h;std::copy(lo+11*i,lo+11*i+11,l.begin());std::copy(hi+11*i,hi+11*i+11,h.begin());
        L* stat=results+4*i;
        stat[3]=split_probe(t,l,h,0,depth,calls,stat)?1:0;
    }
    return 0;
}

// Original certificate preorder remains a complete partition. A zero mask at an
// original leaf is accepted here only if the full new adaptive binary tree passed.
// Previous masks supplied by the caller can encode separately proved exclusions.
struct Adaptive {
    Tester& fan;const std::uint8_t* codes;const std::vector<int>& stop;
    const std::vector<int>& active;std::uint8_t* rem;L* sums;
    int depth_limit,call_limit;
    bool run(int p,Box& lo,Box& hi,int depth){
        if(depth>256)return false;
        int count=active[stop[p]]-active[p];
        if(!count)return true;
        int code=codes[p];
        if(code>=253){
            L result[3]{};
            bool ok=split_probe(fan,lo,hi,0,depth_limit,call_limit,result);
            sums[0]+=result[0];sums[1]=std::max(sums[1],result[1]);
            if(ok){++sums[2];sums[4]+=result[2];}
            else{rem[p]=1;++sums[3];if(result[0]>=call_limit)++sums[5];}
            return true;
        }
        if(count>=64 && depth%5==0){
            ++sums[0];
            if(fan.check(lo,hi)==1){sums[2]+=count;++sums[4];return true;}
        }
        L a=lo[code],b=hi[code];
        if((a+b)%2)return false;
        auto m=(a+b)/2;hi[code]=m;
        if(!run(p+1,lo,hi,depth+1))return false;
        hi[code]=b;lo[code]=m;
        if(!run(stop[p+1],lo,hi,depth+1))return false;
        lo[code]=a;return true;
    }
};
extern "C" int adaptive_trim(const std::uint8_t* codes,int n,const L* lower,const L* upper,
    const L* original_hi,const std::uint8_t* old,std::uint8_t* next,L* stats,
    const std::int32_t* table,const std::int32_t* counts,L goal,int depth,int calls){
    if(n<1 || n>(1<<29) || depth<0 || depth>64 || calls<1 || calls>(1<<24))return 1;
    for(int k=0;k<11;++k){if(lower[k]<=-Fine || upper[k]>=Fine || lower[k]>=upper[k])return 1;}
    std::vector<int> stop(n),active(n+1);
    for(int i=0;i<n;++i){
        if(old[i]&&codes[i]<253)return 2;
        active[i+1]=active[i]+(old[i]?1:0);
    }
    for(int i=n-1;i>=0;--i){
        if(codes[i]>=253)stop[i]=i+1;
        else{if(codes[i]>=11 || i+1>=n || stop[i+1]>=n)return 3;stop[i]=stop[stop[i+1]];}
    }
    if(stop[0]!=n)return 4;
    Tester tester(table,counts,original_hi,goal);
    Adaptive runner{tester,codes,stop,active,next,stats,depth,calls};
    Box lo,hi;std::copy(lower,lower+11,lo.begin());std::copy(upper,upper+11,hi.begin());
    return runner.run(0,lo,hi,0)?0:5;
}
