#include <metal_stdlib>
using namespace metal;

struct Point { float rho,p,a,delta,gamma,c,lambda; uint target,mode; };
struct Row {
    Point p;
    float logFull,logUnit,rate,unitRate,leading,score;
    uint tested,fullPositive,unitPositive,missed,nonfinite,index,salt;
};

inline uint hash32(uint x) {
    x^=x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float unit(uint x,uint k) {
    return float(hash32(x+0x9e3779b9u*k)>>8)*0x1p-24f;
}
inline float logs(uint x,uint k,float a,float b) {return exp2(mix(a,b,unit(x,k)));}

inline float logOneMinusExp(float x) {
    if (x<=0.0f) return -INFINITY;
    if (x>=20.0f) return 0.0f;
    if (x<0.05f)
        return log(x)+log(1.0f+x*(-0.5f+x*(1.0f/6.0f+x*(-1.0f/24.0f+x/120.0f))));
    return log(1.0f-exp(-x));
}

inline void append(thread float2& total,float value) {
    if (value>total.x) {
        total.y=1.0f+total.y*exp(total.x-value);
        total.x=value;
    } else if (isfinite(value)) {
        total.y+=exp(value-total.x);
    }
}
inline float complete(float2 t) { return t.x+log(t.y); }

inline void plateau(thread float2& acc,constant float2* nodes,uint order,
                    float l,float r,float diff,float width,float common,float heat) {
    const float span=r-l;
    if (span<=0.0f) return;
    const float factor=common+logOneMinusExp(heat)+log(span);
    for (uint j=0;j<order;++j) {
        const float u=l+span*nodes[j].x;
        const float dist=diff+width*u;
        append(acc,factor-0.5f*dist*dist+log(nodes[j].y));
    }
}

inline void edge(thread float2& acc,constant float2* nodes,uint order,
                 float near,float far,float sign,float diff,float width,float common,float heat) {
    const float e0=2.0f*(near-0.5f),e1=2.0f*(far-0.5f);
    float length=e1-e0;
    if (length<=0.0f) return;
    const float slope=heat-0.5f*width*(abs(diff)+width);
    if (slope>0) length=min(length,42.0f/slope);
    const float factor=common-heat*e0+log(length*0.5f);
    const float d0=diff+sign*width*near;
    for (uint j=0;j<order;++j) {
        const float z=length*nodes[j].x;
        const float d=d0+sign*0.5f*width*z;
        const float cut=heat*max(0.0f,1.0f-e0-z);
        append(acc,factor-heat*z-0.5f*d*d+logOneMinusExp(cut)+log(nodes[j].y));
    }
}

inline void dip(thread float2& out,constant float2* nodes,uint order,
                float l,float r,float diff,float width,float base,float heat) {
    if (r<=l) return;
    plateau(out,nodes,order,max(l,-0.5f),min(r,0.5f),diff,width,base,heat);
    edge(out,nodes,order,max(0.5f,-r),-l,-1.0f,diff,width,base,heat);
    edge(out,nodes,order,max(0.5f,l),r,1.0f,diff,width,base,heat);
}

kernel void window_scan(constant float2* nodes [[buffer(0)]],device Row* output [[buffer(1)]],
                    constant uint4& params [[buffer(2)]],
                    uint id [[thread_position_in_grid]],uint lane [[thread_index_in_simdgroup]]) {
    const uint index=params.z+id;
    const uint seed=hash32(index^0x51ff3b21u^params.w);
    const uint mode=index%5u,target=2u+hash32(seed+9429u)%5u;
    const float rho=mix(1.75f,4.5f,unit(seed,1));
    const float powers[3]={1.0f,1.5f,2.0f};
    const float p=powers[(index/5u)%3u];
    const float a=logs(seed,2,-7.0f,-2.0f);
    const float delta=logs(seed,3,-8.0f,-1.0f);
    const float gamma=(index/15u)%4u==0u ? -0.5f*delta*unit(seed,4):logs(seed,4,-8.0f,-1.0f);
    float anchor=4.0f;
    for (uint j=0;j<target;++j) anchor*=rho;
    float c=1.0f;
    if (mode==1u) c=mix(1.02f,rho-0.02f,unit(seed,5));
    if (mode==2u) c=sqrt(rho);
    if (mode==3u) c=mix(1.15f,1.35f,unit(seed,5));
    if (mode==4u) c=1.0f+mix(-1.2f,1.2f,unit(seed,5))/anchor;
    const float lambda=c*anchor;
    const float lam2=lambda*lambda,q=0.5f*lam2;
    const Point pt=Point{rho,p,a,delta,gamma,c,lambda,target,mode};
    const float heat=q*(delta+gamma);
    float2 full=float2(-q*delta,1.0f),local=float2(-q*delta-0.3817151463021261f,1.0f);
    float leading=-0.5f*delta,center=4.0f;
    for (uint j=0;j<12u;++j) {
        const float cost=a*powr(center,p);
        const float logw=-2.0794415416798357f-cost;
        const float w=exp(logw),diff=center-lambda;
        const float minDist=max(0.0f,abs(diff)-w);
        const float base=logw+q*gamma-0.9189385332046727f;
        const float upper=base-0.5f*minDist*minDist+0.6931471805599453f;
        leading=max(leading,0.5f*gamma-cost/lam2-0.5f*(diff/lambda)*(diff/lambda));
        if (upper>=full.x-36.0f) {
            dip(full,nodes,params.y,-1.0f,1.0f,diff,w,base,heat);
        }
        if (minDist<=1.0f&&upper>=local.x-36.0f) {
            if (abs(diff)+w<1.0f || (w==0.0f&&abs(diff)<1.0f)) {
                dip(local,nodes,params.y,-1.0f,1.0f,diff,w,base,heat);
            } else if (w==0.0f&&abs(diff)==1.0f) {
                dip(local,nodes,params.y,diff>0.0f?-1.0f:0.0f,
                    diff>0.0f?0.0f:1.0f,diff,w,base,heat);
            } else if (w>0.0f) {
                const float l=max(-1.0f,(-1.0f-diff)/w),r=min(1.0f,(1.0f-diff)/w);
                dip(local,nodes,params.y,l,r,diff,w,base,heat);
            }
        }
        center*=rho;
    }
    const float F=complete(full),U=complete(local),rate=F/lam2,unitRate=U/lam2;
    const bool good=isfinite(F)&&isfinite(U);
    const bool miss=F>0.0f&&U<0.0f;
    const uint g=id/32u;
    for (uint j=0;j<9u;++j) {
        const bool enabled=j==0u || (j==1u&&mode==0u) || (j==2u&&mode==2u)
          || (j==3u&&miss) || (j==4u&&p==2.0f&&lambda>=256.0f)
          || (j==5u&&p<2.0f&&gamma<=0.0f) || (j==6u&&gamma<=0.0f)
          || (j==7u&&lambda>=256.0f) || j==8u;
        float value=abs(rate);
        if (j==3u) value=min(rate,-unitRate);
        if (j==7u) value=abs(rate-leading);
        if (j==8u) value=-abs(F);
        const float val=enabled&&good?value:-INFINITY;
        const float mx=simd_max(val);
        const uint winner=simd_min(val==mx?lane:32u);
        const uint count=simd_sum(uint(enabled));
        const uint pos=simd_sum(uint(enabled&&good&&F>0.0f));
        const uint localPos=simd_sum(uint(enabled&&good&&U>0.0f));
        const uint missCnt=simd_sum(uint(enabled&&good&&miss));
        const uint bad=simd_sum(uint(enabled&&!good));
        if (lane==winner) output[g*9u+j]=Row{pt,F,U,rate,unitRate,leading,mx,count,pos,localPos,missCnt,bad,index,params.w};
    }
}
