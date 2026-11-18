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
// A piecewise linear upper bound on Euclidean length. For t=min/max,
// convexity of sqrt(1+t*t) permits the two secants through the rational
// upper ordinates 1024/1024, 1145/1024, 1449/1024 at t=0,1/2,1.
L norm(L x,L y){
    x=std::abs(x);y=std::abs(y);
    L hi=std::max(x,y),lo=std::min(x,y);
    if(2*lo<=hi)return hi+(242*lo+1023)/1024;
    return (841*hi+608*lo+1023)/1024;
}
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
    std::vector<P> pts,hull,conv,merged;
    std::array<std::vector<int>,5*1025> sort_order;
    static bool before(const P& a,const P& b){
        return a.x!=b.x?a.x<b.x:(a.y!=b.y?a.y<b.y:a.body<b.body);
    }
    L rot_lo[5]{},rot_hi[5]{},width[5]{},hx[5]{},hy[5]{},theta[5]{},tn[5]{},tw[5]{},rad[5]{};
    Tester(const std::int32_t* p,const std::int32_t* n,const L* h,L t):table(p),sizes(n),original_hi(h),goal(t){
        for(int k=0;k<5;++k){if(sizes[k]<3 || sizes[k]>1024)std::abort();total+=sizes[k];}
        if(total>4096)std::abort();
        for(int k=0;k<11;++k)if(h[k]<=-Fine || h[k]>=Fine)std::abort();
        for(int k: {2,5,8})if(h[k]<=0)std::abort();
        merged.resize(total);pts.reserve(total);hull.reserve(2*total);conv.reserve(total);
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
        pts.clear();hull.clear();conv.clear();int off=0,ends[5]{};
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
                    L den=1024*ratio;
                    L lower=1024*lo[k]-anchor,upper=1024*hi[k]-anchor;
                    rot_lo[b]=lower>=0?lower/den:-small(up(-lower,den));
                    rot_hi[b]=upper<=0?-((-upper)/den):small(up(upper,den));
                    width[b]=small(up(hi[k]-lo[k],2*ratio));
                    theta[b]=std::max(std::abs(rot_lo[b]),std::abs(rot_hi[b]));
                }
            }
            tw[b]=norm(hx[b],hy[b]);tn[b]=norm(tx,ty);
            const auto* block=table+2*(row*total+off);
            auto& order=sort_order[1025*b+row];
            if(order.empty()){
                order.resize(sizes[b]);
                for(int j=0;j<sizes[b];++j) order[j]=j;
                std::sort(order.begin(),order.end(),[&](int u,int v){
                    auto x=block+2*u,y=block+2*v;
                    return x[0]!=y[0]?x[0]<y[0]:x[1]<y[1];
                });
            }
            for(int j:order){
                L a=block[2*j],c=block[2*j+1];
                L x=tx+a,y=ty+c;
                // Keeping any subset of valid witnesses is always sound. This
                // threshold removes points well inside the fixed disk-sample hull.
                constexpr L cut=Q*499/1000;
                if(b && x*x+y*y<cut*cut)continue;
                pts.push_back(P{x,y,a,c,b});
            }
            off+=sizes[b];ends[b]=int(pts.size());
        }
        // Translation preserves each body's lexicographic ordering. Merge its
        // nominal angle row with the other four exactly sorted lists.
        merged.resize(pts.size());auto a=pts.begin(),b=merged.begin();
        int j1=ends[0],j2=ends[1],j3=ends[2],j4=ends[3],npts=ends[4];
        std::merge(a,a+j1,a+j1,a+j2,b,before);
        std::merge(a+j2,a+j3,a+j3,a+j4,b+j2,before);
        
        std::merge(b,b+j2,b+j2,b+j4,a,before);
        std::merge(a,a+j4,a+j4,a+npts,b,before);
        pts.swap(merged);
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
        W twice=0,error=0,grot[5]{},cosine[5]{},cross_sum[5][5]{},dot_sum[5][5]{};
        V g[5]{},mixed[5][5]{};
        int translation_pairs[5][5]{};
        auto n=conv.size();
        for(unsigned i=0;i<n;++i){
            const P& a=conv[i];const P& next=conv[(i+1)%n];const P& prev=conv[(i+n-1)%n];
            L dx=next.x-prev.x,dy=next.y-prev.y;
            twice+=det(a.x,a.y,next.x,next.y);
            twice+=curved_gain(a,next);
            int k=a.body;g[k].x+=dx;g[k].y+=dy;
            grot[k]-=W(a.a)*dx+W(a.b)*dy;
            int j=next.body;
            cosine[k]+=det(a.a,a.b,dx,dy);
            mixed[k][j].x+=next.a;mixed[k][j].y+=next.b;
            mixed[j][k].x-=a.a;mixed[j][k].y-=a.b;
            W cross=det(a.a,a.b,next.a,next.b);
            if(k==j){
                // det((R-I)a,(R-I)b)=2(1-cos(alpha))*det(a,b).
                cosine[k]-=2*cross;
            }else{
                int first=std::min(k,j),last=std::max(k,j),sign=k<j?1:-1;
                translation_pairs[first][last]+=sign;
                cross_sum[first][last]+=cross;
                dot_sum[first][last]+=sign*(W(a.a)*next.a+W(a.b)*next.b);
            }
        }
        for(int k=1;k<5;++k){
            error+=W(std::abs(g[k].y))*hx[k]+W(std::abs(g[k].x))*hy[k];
            // Actual interval displacement from its witness row can be one-sided.
            // |sin(a)-a| <= |a|^3/6 for every real a.
            W scalar=grot[k]*(grot[k]>=0?rot_lo[k]:rot_hi[k]);
            twice+=scalar>=0?scalar/Q:-up(-scalar,Q);
            W first=up(ab(grot[k])*theta[k],Q);
            error+=std::min(first,up(up(first*theta[k],Q)*theta[k],6*Q));
            if(cosine[k]>0)error+=up(up(cosine[k]*theta[k],Q)*theta[k],2*Q);
            for(int j=2;j<5;++j){
                auto v=mixed[k][j];
                if(!v.x&&!v.y)continue;
                W cs=W(hx[k])*std::abs(v.y)+W(hy[k])*std::abs(v.x);
                W sn=W(hx[k])*std::abs(v.x)+W(hy[k])*std::abs(v.y);
                W component=up(sn*theta[j],Q)+up(up(cs*theta[j],Q)*theta[j],2*Q);
                W radial=W(tw[k])*up(W(norm(v.x,v.y))*theta[j],Q);
                error+=std::min(component,radial);
            }
            for(int j=k+1;j<5;++j){
                W box_det=W(hx[k])*hy[j]+W(hy[k])*hx[j];
                error+=std::abs(translation_pairs[k][j])*std::min(box_det,W(tw[k])*tw[j]);
                if(!theta[k] || !theta[j])continue;
                W cross=up(up(ab(cross_sum[k][j])*theta[k],Q)*theta[j],Q);
                W dot=up(up(ab(dot_sum[k][j])*theta[k],Q)*theta[j],Q);
                error+=cross+std::min(dot,up(dot*(theta[k]+theta[j]),2*Q));
            }
        }
        W best=-1;
        for(int b=1;b<5;++b){
            int k=3*(b-1);
            auto consider=[&](int index,W score){if(hi[index]-lo[index]>=2 && score>best){split_axis=index;best=score;}};
            consider(k,W(hx[b])*(std::abs(g[b].y)+Q/30));
            consider(k+1,W(hy[b])*(std::abs(g[b].x)+Q/30));
            if(b>1)consider(k-1,W(width[b])*(ab(grot[b])/Q+Q/100));
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
    int k=t.split_axis;auto l=lo[k],h=hi[k];if(h-l<2)return false;
    auto mid=l+(h-l)/2;hi[k]=mid;
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
    if(n<1 || n>(1<<29) || depth<0 || depth>128 || calls<1 || calls>(1<<24))return 1;
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
