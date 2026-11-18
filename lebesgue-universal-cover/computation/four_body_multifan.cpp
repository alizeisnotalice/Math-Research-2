// Integer fan-area lower bound on whole 11D dyadic boxes; coordinates in units 2^-28.
#include <algorithm>
#include <array>
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
struct Form { W f=0; std::array<W,8> g{}; };
using Move=std::array<int,8>;
struct Tester {
    const std::int32_t* table;const std::int32_t* sizes;const L* original_hi;
    L goal;int total=0;int split_axis=0;int modes=4;
    std::vector<P> pts,hull,conv;
    L hx[5]{},hy[5]{},theta[5]{},tn[5]{},tw[5]{},rad[5]{};
    Tester(const std::int32_t* p,const std::int32_t* n,const L* h,L t):table(p),sizes(n),original_hi(h),goal(t){
        for(int k=0;k<5;++k){if((k<4 && sizes[k]<3) || (k==4 && sizes[k]!=0) || sizes[k]>1024)std::abort();total+=sizes[k];}
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
    int form(const Box& lo,const Box& hi,const Move& move,Form& result){
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
        // Integer guide positions choose material and ordering only. Validation below
        // uses the original exact anchors on the full input box.
        const bool shifted = std::any_of(move.begin(),move.end(),[](int v){return v!=0;});
        auto guide = [&](const P& p) -> P {
            if(!shifted||!p.body)return p;
            const int i=3*(p.body-1);
            const L tx=move[i]*hx[p.body], ty=move[i+1]*hy[p.body];
            const L da=i?move[i-1]*theta[p.body]:0;
            return P{p.x+tx-L(W(p.b)*da/Q),p.y+ty+L(W(p.a)*da/Q),p.a,p.b,p.body};
        };
        if(shifted)for(auto& p:pts)p=guide(p);
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
        if(shifted)for(auto& p:conv)if(p.body){
            const int i=3*(p.body-1);
            const L da=i?move[i-1]*theta[p.body]:0;
            p.x-=move[i]*hx[p.body]-L(W(p.b)*da/Q);
            p.y-=move[i+1]*hy[p.body]+L(W(p.a)*da/Q);
        }
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
        int turns=0;
        for(unsigned j=0;j<conv.size();++j)
            if(conv[j].y<=0 && conv[(j+1)%conv.size()].y>0)++turns;
        if(turns!=1)return 0;
        W twice=0,error=0,grot[5]{};V g[5]{};
        auto n=conv.size();
        for(unsigned i=0;i<n;++i){
            const P& a=conv[i];const P& next=conv[(i+1)%n];const P& prev=conv[(i+n-1)%n];
            L dx=next.x-prev.x,dy=next.y-prev.y;
            twice+=det(a.x,a.y,next.x,next.y);
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
        result.f=twice-error;
        for(int k=1;k<4;++k){
            int i=(k-1)*3;
            result.g[i]=W(g[k].y)*hx[k];
            result.g[i+1]=-W(g[k].x)*hy[k];
            if(i){W v=grot[k]*theta[k]; result.g[i-1]=v/Q;
                if(v%Q) --result.f;
            }
        }
        W best=-1;
        for(int b=1;b<4;++b){
            int k=3*(b-1);
            auto consider=[&](int index,W score){if(score>best){split_axis=index;best=score;}};
            consider(k,W(hx[b])*(std::abs(g[b].y)+Q/30));
            consider(k+1,W(hy[b])*(std::abs(g[b].x)+Q/30));
            if(b>1)consider(k-1,W(theta[b])*(ab(grot[b])/Q+Q/100));
        }
        return 1;
    }
    static W cost(const Form& a){
        W z=a.f;for(W g:a.g)z-=ab(g);return z;
    }
    static Form blend(const Form& a,const Form& b,L w){
        // Rational mixture rounded conservatively to the same integer denominator.
        constexpr L den=1024;Form out;W z=(den-w)*a.f+w*b.f;
        out.f=z/den-(z%den<0);
        for(int j=0;j<8;++j){
            W t=(den-w)*a.g[j]+w*b.g[j];
            out.g[j]=t/den;if(t%den)--out.f;
        }return out;
    }
    Form bestpair(const Form& a,const Form& b){
        Form best=cost(a)>cost(b)?a:b;
        for(int j=-1;j<8;++j){
            L w=512;
            if(j>=0){const W ap=a.g[j],bp=b.g[j];
                if((ap>=0&&bp>=0)||(ap<=0&&bp<=0))continue;
                w=L((1024*ap+(ap-bp)/2)/(ap-bp));
            }
            if(w<=0||w>=1024)continue;
            const Form c=blend(a,b,w);if(cost(c)>cost(best))best=c;
        }return best;
    }
    int check(const Box& lo,const Box& hi){
        Move move{};Form f;
        if(form(lo,hi,move,f)!=1)return 0;
        if(cost(f)>=goal)return 1;
        const int original_axis=split_axis;
        std::array<Form,12> forms{};forms[0]=f;int used=1;
        Form best=f;
        for(int i=0;i<modes;++i){
            move={};
            for(int j=0;j<8;++j){
                int sign=best.g[j]<0?1:best.g[j]>0?-1:0;
                if(i==0)move[j]=sign;
                if(i==1)move[j]=f.g[j]<0?-1:f.g[j]>0?1:0;
                if(i==2)move[j]=sign;
                if(i==3 && j==original_axis)move[j]=sign;
                if(i>=4 && j==(original_axis+(i-3)/2)%8)move[j]=i%2?1:-1;
            }
            Form next;
            if(form(lo,hi,move,next)!=1)continue;
            if(cost(next)>=goal)return 1;
            Form improved=bestpair(best,next);
            for(int k=0;k<used;++k){
                Form c=bestpair(forms[k],next);if(cost(c)>cost(improved))improved=c;
            }
            if(cost(improved)>=goal)return 1;
            best=improved;forms[used++]=next;
        }
        // Restore the original fan's proven split heuristic after guide generation.
        split_axis=original_axis;return 0;
    }
};
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

