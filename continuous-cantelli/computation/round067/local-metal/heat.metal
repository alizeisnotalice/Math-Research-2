#include <metal_stdlib>
using namespace metal;

struct Point {
    float c1,s1,c2,s2,c3,s3, b1,t1,b2,t2,b3,t3, d,eta,p;
    uint g,n;
};
struct Row {
    Point p;
    float real,imag,linearReal,linearImag,reference,score;
    uint tested,nonfinite,index,salt;
};

inline uint hash32(uint x) {
    x ^= x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float unit(uint x,uint k) {
    return float(hash32(x+0x9e3779b9u*k)>>8)*0x1p-24f;
}
inline float expminus(float x) {
    return x*(1.0f+x*(0.5f+x*(1.0f/6.0f+x*(1.0f/24.0f+x*(1.0f/120.0f+x/720.0f)))));
}

inline float harmonic(float x,float c1,float s1,float c2,float s2,float c3,float s3) {
    const float t=cos(x), u=sin(x);
    const float t2=t*t-u*u, u2=2.0f*t*u;
    return c1*t+s1*u+c2*t2+s2*u2+c3*(t2*t-u2*u)+s3*(u2*t+t2*u);
}

inline float4 panel(Point p,constant float2* nodes,uint order,
                    float lo,float hi,float origin,float phase,bool scaled) {
    const float width=hi-lo;
    if (width<=0) return 0;
    const bool first=lo==0.0f;
    const float gf=float(p.g), nf=float(p.n);
    float4 accum=0;
    for (uint j=0;j<order;++j) {
        const float z=nodes[j].x;
        const float t=lo+width*(first ? z*z*z*z : z);
        const float w=width*nodes[j].y*(first ? 4.0f*z*z*z : 1.0f);
        const float radial=2.0f*sin(0.5f*gf*t)/gf;
        const float rate=0.5f*p.d*nf*nf*powr(radial,p.p);
        for (int side=-1;side<=1;side+=2) {
            const float u=float(side)*t;
            const float q0=1.0f+harmonic(gf*u,p.c1,p.s1,p.c2,p.s2,p.c3,p.s3);
            const float qb=harmonic(origin+u,p.b1,p.t1,p.b2,p.t2,p.b3,p.t3);
            const float base=exp(-rate*q0);
            float primary,linear=0;
            if (scaled) {
                const float shift=-rate*p.eta*qb;
                if (abs(shift)<0.1f) {
                    primary=base*expminus(shift)/p.eta;
                } else {
                    primary=(exp(-rate*(q0+p.eta*qb))-base)/p.eta;
                }
                linear=-rate*qb*base;
            } else {
                primary=exp(-rate*(q0+p.eta*qb));
            }
            const float angle=phase+nf*u;
            const float2 direction=float2(cos(angle),sin(angle));
            accum+=w*float4(primary*direction,linear*direction);
        }
    }
    return accum;
}

inline float4 coefficient(Point p,constant float2* nodes,uint order) {
    const bool off=p.n%p.g!=0;
    if (off && p.eta==0.0f) return 0;
    const float pi=3.14159265358979323846f;
    const float edge=pi/float(p.g);
    const float width=powr(2.0f/(p.d*float(p.n*p.n)),1.0f/p.p);
    const float end=min(edge,0.5f*pi*width*powr(37.0f/0.35f,1.0f/p.p));
    const float cuts[5]={0,min(end,width),min(end,width*powr(4.0f,1.0f/p.p)),
                         min(end,width*powr(16.0f,1.0f/p.p)),end};
    float4 total=0;
    for (uint g=0;g<p.g;++g) {
        const float origin=(2.0f*pi/float(p.g))*float(g);
        const float phase=(2.0f*pi/float(p.g))*float((p.n*g)%p.g);
        for (uint j=0;j<4;++j) {
            const float span=cuts[j+1]-cuts[j];
            const uint count=max(1u,uint(ceil(span*float(p.n)/22.0f)));
            const float step=span/float(count);
            for (uint k=0;k<count;++k) {
                total+=panel(p,nodes,order,cuts[j]+float(k)*step,
                             cuts[j]+float(k+1)*step,origin,phase,off);
            }
        }
    }
    return total*(0.5f/pi);
}

kernel void heat(constant float2* nodes [[buffer(0)]],device Row* output [[buffer(1)]],
                 constant uint4& params [[buffer(2)]],
                 uint id [[thread_position_in_grid]],uint lane [[thread_index_in_simdgroup]]) {
    const uint i=params.z+id;
    const uint r=hash32(i^0x0b6aa1fbu^params.w);
    const uint pg=i%11u;
    const float powers[11]={0.5f,0.75f,1,1.5f,1.75f,2,2.25f,3,4,6,8};
    const float gamma[5]={2.0f,1.19063934876f,1.0f,0.902745292951f,0.890617733087f};
    float a[12], sum0=0,sum1=0;
    for (uint k=0;k<12;++k) {
        a[k]=2.0f*unit(r,k+1)-1.0f;
        if (k<6) sum0+=abs(a[k]); else sum1+=abs(a[k]);
    }
    for (uint k=0;k<12;++k) a[k]*=(k<6?0.4f/(sum0+0.2f):1.0f/(sum1+0.2f));
    const uint g=1u+(i/11u)%4u;
    const uint n=1u+hash32(r+942771u)%96u;
    const float d=exp2(mix(-6.64385619f,4.0f,unit(r,13)));
    const float eta=(i/44u)%3u==0u ? 0.0f : exp2(mix(-18.0f,-2.0f,unit(r,14)));
    const Point p=Point{a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],
                        d,eta,powers[pg],g,n};
    const float4 values=coefficient(p,nodes,params.y);
    const bool off=n%g!=0, single=g==1u;
    const float q0=1.0f+p.c1+p.c2+p.c3+eta*(p.b1+p.b2+p.b3);
    const float rate=d*q0;
    float reference=0;
    if (single && pg<5u) {
        reference=gamma[pg]*0.3183098861837907f*powr(2.0f/(rate*float(n*n)),1.0f/p.p);
    }
    if (single && pg==5u) {
        reference=exp(-0.5f/rate)*rsqrt(6.283185307179586f*rate)/float(n);
    }
    const bool good=all(isfinite(values));
    const float primary=length(values.xy);
    const uint group=id/32u;
    for (uint j=0;j<10u;++j) {
        const bool enabled=j==0u || (j==1u && single && pg<5u && n>=24u && rate>=0.25f)
            || (j==2u && single && pg==5u && n>=24u && rate>=0.25f)
            || (j==3u && single && pg>5u && n>=24u)
            || (j==4u && off && eta==0.0f)
            || (j==5u && off && eta>0.0f)
            || (j==6u && off && eta>0.0f && eta<=0.001f)
            || (j==7u && !single && !off)
            || (j==8u && pg<=2u && n>=64u && single)
            || (j==9u && pg>=8u && n>=64u && single);
        float metric=primary;
        if ((j==1u||j==2u) && enabled) metric=length(values.xy-float2(reference,0.0f))/reference;
        if (j==6u && enabled) metric=length(values.xy-values.zw);
        const float val=enabled&&good?metric:-INFINITY;
        const float mx=simd_max(val);
        const uint winner=simd_min(val==mx?lane:32u);
        const uint count=simd_sum(uint(enabled));
        const uint bad=simd_sum(uint(enabled&&!good));
        if (lane==winner) output[group*10u+j]=Row{p,values.x,values.y,values.z,values.w,reference,mx,count,bad,i,params.w};
    }
}
