#include <metal_stdlib>
using namespace metal;

struct Point { float beta,r,g,amp,t; uint family; };
struct Row {
    Point p;
    float positive,negative,logJPositive,logJNegative,logMeanPositive,logMeanNegative,score;
    uint tested,positiveAbove,negativeAbove,opposite,nonfinite,index,salt;
};

inline uint hash32(uint x) {
    x ^= x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float unit(uint x,uint k) {
    return float(hash32(x+0x9e3779b9u*k)>>8)*0x1p-24f;
}
inline float log_phi2(float z) {
    const float a=abs(z);
    if (a<=0.3f)
        return log(0.5f+z*(1.0f/6.0f+z*(1.0f/24.0f+z*(1.0f/120.0f+z*(1.0f/720.0f+
                   z*(1.0f/5040.0f+z/40320.0f))))));
    const float normal=2.0f*log(a);
    if (z>=12.0f) return z+log(1.0f-(1.0f+z)*exp(-z))-normal;
    return log(exp(z)-1.0f-z)-normal;
}
inline float logmean(Point p,float g) {
    const float B=1.0f+p.beta;
    const float xi=p.t*rsqrt(B);
    return log(2.0f*p.amp)-1.5f*log(B)-0.5f*p.beta/B*p.t*p.t
          +log(1.0f+g*xi+p.r*xi*xi);
}
inline void update(thread float2& acc,float v) {
    if (v>acc.x) {
        acc.y=1.0f+acc.y*exp(acc.x-v);
        acc.x=v;
    } else if (isfinite(v)) {
        acc.y+=exp(v-acc.x);
    }
}
inline void panel(Point p,constant float2* nodes,uint order,float lo,float hi,
                  thread float2& plus,thread float2& minus) {
    const float width=hi-lo;
    if (width<=0) return;
    const float root=sqrt(1.0f+p.beta);
    const float s=0.5f*p.t*p.t;
    const float logamp=log(p.amp);
    for (uint j=0;j<order;++j) {
        const float x=lo+width*nodes[j].x;
        const float y=root*x,y2=y*y;
        const float even=1.0f-y2-p.r*(y2*y2-6.0f*y2+3.0f);
        const float odd=-p.g*y*(y2-3.0f);
        const float exponent=logamp-0.5f*p.beta*x*x;
        const float d=x-p.t;
        const float factor=log(width*nodes[j].y)-0.5f*d*d-0.9189385332046727f;
        const float P=even+odd;
        if (P!=0.0f) {
            const float b=exponent+log(abs(P));
            const float z=copysign(s*exp(b),P);
            update(plus,factor+2.0f*b+log_phi2(z));
        }
        if (p.g!=0.0f) {
            const float Pm=even-odd;
            if (Pm!=0.0f) {
                const float b=exponent+log(abs(Pm));
                const float z=copysign(s*exp(b),Pm);
                update(minus,factor+2.0f*b+log_phi2(z));
            }
        }
    }
}

inline float2 evaluate(Point p,constant float2* nodes,uint order) {
    const float w=rsqrt(1.0f+p.beta), tiny=w/max(1.0f,p.t);
    float cuts[15]={p.t-24.0f,p.t+24.0f,p.t-6.0f,p.t-2.0f,p.t,p.t+2.0f,p.t+6.0f,
                    -5.0f*w,-w,-tiny,0,tiny,w,5.0f*w,12.0f*w};
    for (uint j=1;j<15;++j) {
        float x=cuts[j];
        uint k=j;
        while (k>0&&x<cuts[k-1]) { cuts[k]=cuts[k-1]; --k; }
        cuts[k]=x;
    }
    float2 plus=float2(-INFINITY,0.0f),minus=float2(-INFINITY,0.0f);
    for (uint j=0;j+1<15;++j) {
        panel(p,nodes,order,max(cuts[j],p.t-24.0f),min(cuts[j+1],p.t+24.0f),plus,minus);
    }
    const float P=plus.x+log(plus.y);
    return float2(P,p.g==0.0f?P:minus.x+log(minus.y));
}

kernel void jensen(constant float2* nodes [[buffer(0)]],device Row* output [[buffer(1)]],
                   constant uint4& params [[buffer(2)]],
                   uint id [[thread_position_in_grid]],uint lane [[thread_index_in_simdgroup]]) {
    const uint index=params.z+id;
    const uint r=hash32(index^0x9b2f31ebu^params.w);
    const uint family=index%3u;
    const float beta=exp2(mix(-6.0f,9.0f,unit(r,1)));
    const float h=family==0u ? 0.0f : unit(r,2)/12.0f;
    const float g=family<2u ? 0.0f : (unit(r,3)-0.5f)*sqrt(h);
    const float A=exp2(mix(-10.0f,-0.02914635f,unit(r,4)));
    const float t=(index/3u)%3u==0u ? mix(4.0f,12.0f,unit(r,5)) : exp2(mix(-10.0f,3.5849625f,unit(r,5)));
    const Point p=Point{beta,h,g,A,t,family};
    const float2 j=evaluate(p,nodes,params.y);
    const float m=logmean(p,g),mm=logmean(p,-g);
    const float positive=j.x-m,negative=j.y-mm;
    const bool good=isfinite(positive)&&isfinite(negative);
    const bool opposite=(positive>0)!=(negative>0);
    const bool even=g==0.0f;
    const float threshold=A*(1.0f-3.0f*h)*(1.0f+beta);
    const uint group=id/32u;
    for (uint s=0;s<9;++s) {
        const bool enabled=s==0 || (s==1&&t<=0.05f)
          || (s==2&&even&&t>=6.0f&&threshold<0.95f)
          || (s==3&&even&&t>=6.0f&&threshold>1.05f)
          || (s==4&&!even&&t>=6.0f)
          || (s==5&&!even&&opposite)
          || s==6 || (s==7&&A<0.1f) || (s==8&&even&&t>=6.0f&&abs(threshold-1.0f)<=0.05f);
        float score=max(abs(positive),abs(negative));
        if (s==5) score=min(abs(positive),abs(negative));
        if (s==6) score=-min(abs(positive),abs(negative));
        if (s==2) score=-positive;
        if (s==3) score=positive;
        const float val=enabled&&good?score:-INFINITY;
        const float best=simd_max(val);
        const uint winner=simd_min(best==val?lane:32u);
        const uint count=simd_sum(uint(enabled));
        const uint above1=simd_sum(uint(enabled&&good&&positive>0.0f));
        const uint above2=simd_sum(uint(enabled&&good&&negative>0.0f));
        const uint split=simd_sum(uint(enabled&&good&&opposite));
        const uint bad=simd_sum(uint(enabled&&!good));
        if (lane==winner)
            output[group*9u+s]=Row{p,positive,negative,j.x,j.y,m,mm,best,count,above1,above2,split,bad,index,params.w};
    }
}
