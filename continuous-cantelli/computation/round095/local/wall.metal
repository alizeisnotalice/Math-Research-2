#include <metal_stdlib>
using namespace metal;

struct Point { float k,v; uint family; };
struct Row {
    Point p;
    float logR,logLinear,oldDelta,upperDelta,lowerDelta,score;
    uint tested,oldAbove,upperAbove,nonfinite,index,salt;
};

inline uint hash32(uint x) {
    x^=x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float unit(uint x,uint k) {
    return float(hash32(x+0x9e3779b9u*k)>>8)*0x1p-24f;
}
inline void update(thread float2 &a,float y) {
    if (y>a.x) {
        a.y=1+a.y*exp(a.x-y);
        a.x=y;
    } else if (isfinite(y)) {
        a.y+=exp(y-a.x);
    }
}

inline float exponent(float k,float v,float w) {
    const float z=v*w;
    const float a=k-v,center=max(-a,0.0f);
    if (z>=12.0f) {
        const float d=w-center;
        return -0.5f*d*d-max(a,0.0f)*w-log(2.0f*z)+log(1.0f-exp(-2.0f*z));
    }
    float logsinc;
    if (z<=0.25f) {
        const float x=z*z;
        logsinc=log(1.0f+x*(1.0f/6.0f+x*(1.0f/120.0f+x*(1.0f/5040.0f+x/362880.0f))));
    } else {
        logsinc=log(sinh(z)/z);
    }
    return -0.5f*(w*w+center*center)-k*w+logsinc;
}

inline void panel(float k,float v,constant float2 *nodes,uint n,float a,float b,
                  thread float2 &total) {
    const float width=b-a;
    if (width<=0) return;
    for (uint j=0;j<n;++j) {
        const float w=a+width*nodes[j].x;
        update(total,log(width*nodes[j].y*w)+exponent(k,v,w));
    }
}

inline float integrate(float k,float v,constant float2 *nodes,uint n) {
    const float a=k-v,center=max(0.0f,-a);
    float cuts[9];
    if (a>=0.0f) {
        const float scale=1.0f/max(1.0f,a),end=min(14.0f,48.0f/max(a,0.01f));
        cuts[0]=0; cuts[1]=min(end,scale); cuts[2]=min(end,4.0f*scale);
        cuts[3]=min(end,12.0f*scale); cuts[4]=min(end,24.0f*scale);
        cuts[5]=end; cuts[6]=end; cuts[7]=end; cuts[8]=end;
    } else {
        cuts[0]=max(0.0f,center-14.0f); cuts[1]=max(0.0f,center-8.0f);
        cuts[2]=max(0.0f,center-4.0f); cuts[3]=max(0.0f,center-1.0f);
        cuts[4]=center; cuts[5]=center+1.0f; cuts[6]=center+4.0f;
        cuts[7]=center+8.0f; cuts[8]=center+14.0f;
    }
    float2 total=float2(-INFINITY,0.0f);
    for (uint j=0;j<8;++j) panel(k,v,nodes,n,cuts[j],cuts[j+1],total);
    return total.x+log(total.y)+0.5f*center*center;
}

kernel void wall(constant float2* nodes [[buffer(0)]],device Row* output [[buffer(1)]],
                 constant uint4& params [[buffer(2)]],
                 uint id [[thread_position_in_grid]],uint lane [[thread_index_in_simdgroup]]) {
    const uint index=params.z+id;
    const uint seed=hash32(index^0x28f19025u^params.w);
    const uint family=index%6u;
    float k=exp2(mix(-12.0f,8.0f,unit(seed,1)));
    float v=exp2(mix(-24.0f,6.0f,unit(seed,2)));
    if (family==0u) v=0;
    if (family==2u) v=min(64.0f,k*mix(0.7f,1.2f,unit(seed,3)));
    if (family==3u) {
        v=exp2(mix(-3.0f,6.0f,unit(seed,4)));
        k=max(0x1p-12f,v+mix(-8.0f,8.0f,unit(seed,5)));
    }
    if (family==4u) {
        k=exp2(mix(2.0f,8.0f,unit(seed,6)));
        v=0.1f*k*unit(seed,7);
    }
    if (family==5u) {
        v=exp2(mix(0.0f,6.0f,unit(seed,8)));
        k=min(k,v);
    }
    const Point p=Point{k,v,family};
    const float a=k-v;
    const float r=integrate(k,v,nodes,params.y);
    const float lin=v==0.0f?r:integrate(k,0,nodes,params.y);
    const float old=r+2.0f*log(k);
    const float ceiling=a>0.0f?r+log(a)+log(k+v):0.0f;
    const float lower=r-lin;
    const bool good=isfinite(r)&&isfinite(lin)&&isfinite(old)&&isfinite(ceiling);
    const uint group=id/32u;
    for (uint j=0;j<9;++j) {
        const bool enabled=j==0u || (j==1u&&k>=16.0f&&v/k<=0.1f)
            || (j==2u&&a>=4.0f) || (j==3u&&a<=-4.0f)
            || (j==4u&&abs(a)<=0.125f&&k>=0.5f) || (j==5u&&v<=1e-5f)
            || (j==6u&&a>0.0f&&old>0.0f) || (j==7u&&a>0.0f&&v/k>=0.8f)
            || (j==8u&&k<=0.25f&&v>=1.0f);
        float metric=abs(r);
        if (j==1u) metric=abs(old);
        if (j==2u) metric=ceiling;
        if (j==5u) metric=abs(lower);
        if (j==6u||j==7u) metric=old;
        const float val=enabled&&good?metric:-INFINITY;
        const float best=simd_max(val);
        const uint winner=simd_min(val==best?lane:32u);
        const uint count=simd_sum(uint(enabled));
        const uint oldCount=simd_sum(uint(enabled&&good&&old>0.0f));
        const uint upperCount=simd_sum(uint(enabled&&good&&a>0.0f&&ceiling>0.0f));
        const uint bad=simd_sum(uint(enabled&&!good));
        if (lane==winner)
            output[group*9u+j]=Row{p,r,lin,old,ceiling,lower,best,count,oldCount,upperCount,bad,index,params.w};
    }
}
