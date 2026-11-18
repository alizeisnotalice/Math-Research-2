#include <metal_stdlib>
using namespace metal;
constant float tau=6.2831853071795864769f;
constant uint powers[15]={1u,3u,9u,27u,81u,243u,729u,2187u,6561u,19683u,59049u,177147u,531441u,1594323u,4782969u};

inline int position(device const uint* packed,uint n,uint global,uint ix) {
    const uint scale=powers[global-n],at=(1u<<n)-1u+ix;
    return int(packed[at]*scale+(scale-1u)/2u);
}
inline uint lower(device const uint* pts,uint n,uint N,int target) {
    uint a=0,b=1u<<n;
    while(a<b) {
        const uint m=(a+b)/2u;
        if(position(pts,n,N,m)<target) a=m+1u;
        else b=m;
    }
    return a;
}
inline float2 term(int q,int x,float offset,float unit,float invSigma,int period,float ga) {
    const float t=(float(q-x)*unit-offset)*invSigma;
    const float wing=2.0f+(float(x)+.5f)*unit+offset;
    const float weight=ga*exp(-.5f*((wing-2.0f)*(wing+2.0f)+t*t));
    int frac=(q-x)%period;
    if(frac>period/2) frac-=period;
    if(frac<-(period/2)) frac+=period;
    const float theta=tau*(float(frac)/float(period)-offset/(float(period)*unit));
    return weight*float2(cos(theta),-sin(theta));
}
kernel void cantor(device const uint4* tasks [[buffer(0)]],
                   device const float* sigmas [[buffer(1)]],
                   device const uint* positions [[buffer(2)]],
                   constant float2* nw [[buffer(3)]],
                   constant uint& order [[buffer(4)]],
                   device float* result [[buffer(5)]],
                   uint id [[thread_position_in_grid]]) {
    const uint4 item=tasks[id];
    const uint n=item.x,j=item.y,m=item.z,N=max(n,j);
    const float unit=1.0f/float(powers[N]),sig=sigmas[id],invSigma=1.0f/sig;
    const int a=position(positions,n,N,item.w),period=int(powers[N-m]);
    const float pa=2.0f+(float(a)+.5f)*unit;
    const float ga=exp(-.5f*(pa-2.0f)*(pa+2.0f));
    const int radius=int(ceil(10.0f*sig/unit));
    float2 E=0.0f,H=0.0f,G=0.0f;
    const uint l=lower(positions,n,N,a-radius),r=lower(positions,n,N,a+radius+1);
    for(uint i=l;i<r;++i) {
        const int x=position(positions,n,N,i);
        E+=term(a,x,0.0f,unit,invSigma,period,ga);
    }
    E*=1.0f/float(1u<<n);
    const float gapLen=1.0f/float(powers[j+1u]);
    const int spread=int(ceil(0.5f*gapLen/unit));
    const uint gl=lower(positions,j,N,a-radius-spread);
    const uint gr=lower(positions,j,N,a+radius+spread+1);
    for(uint i=gl;i<gr;++i) {
        const int x=position(positions,j,N,i);
        if(abs(a-x)<=radius) H+=term(a,x,0.0f,unit,invSigma,period,ga);
        for(uint k=0;k<order;++k)
            G+=nw[k].y*term(a,x,nw[k].x*gapLen,unit,invSigma,period,ga);
    }
    H*=1.0f/float(1u<<j);
    G*=1.0f/float(1u<<j);
    for(uint k=0;k<2u;++k) {
        result[id*6u+k]=E[k];
        result[id*6u+2u+k]=H[k];
        result[id*6u+4u+k]=G[k];
    }
}
