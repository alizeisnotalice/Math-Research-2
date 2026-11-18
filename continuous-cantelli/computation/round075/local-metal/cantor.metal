#include <metal_stdlib>
using namespace metal;

struct Point { float lambda,p,d; uint n,ternary,frequency; };
struct Result { float full,coarse,fine,tail,approximation; uint levels; };
struct Row { Point p; float full,coarse,fine,tail,approximation,score; uint levels,tested,nonfinite,index,salt; };

inline uint hash32(uint x) {
    x ^= x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float unit(uint x,uint k) {
    return float(hash32(x+0x9e3779b9u*k)>>8)*0x1p-24f;
}
inline float expminus(float x) {
    if (x < -0.06f) return exp(x)-1.0f;
    return x*(1.0f+x*(0.5f+x*(1.0f/6.0f+x*(1.0f/24.0f+x*(1.0f/120.0f+x/720.0f)))));
}

inline float gap_panel(constant float2* nodes,uint order,float lo,float hi,float phase,
                       float omega,float s,float p,bool subtract) {
    const float width=hi-lo;
    if (width<=0) return 0;
    const bool first=lo==0.0f;
    float total=0;
    for (uint k=0;k<order;++k) {
        const float z=nodes[k].x,z2=z*z;
        const float t=lo+width*(first ? z2*z2 : z);
        const float jac=width*(first ? 4.0f*z2*z : 1.0f);
        const float exponent=-s*powr(t,p);
        total+=nodes[k].y*jac*cos(phase-omega*t)*(subtract ? expminus(exponent):exp(exponent));
    }
    return total;
}

inline float gap_integral(constant float2* nodes,uint order,float phase,float omega,float s,float p) {
    const bool subtract=s<=0.5f;
    const float scale=powr(1.0f/s,1.0f/p);
    const float cap=subtract?1.0f:min(1.0f,scale*powr(42.0f,1.0f/p));
    const float cuts[4]={0,min(cap,scale),min(cap,scale*powr(8.0f,1.0f/p)),cap};
    float total=0.0f;
    for (uint k=0;k<3;++k) {
        const float w=cuts[k+1]-cuts[k];
        if (w<=0) continue;
        const uint count=max(1u,uint(ceil(w*omega/20.0f)));
        const float step=w/float(count);
        for (uint j=0;j<count;++j) {
            total+=gap_panel(nodes,order,cuts[k]+float(j)*step,cuts[k]+float(j+1)*step,
                            phase,omega,s,p,subtract);
        }
    }
    if (!subtract) {
        float sinc;
        if (omega<0.01f) {
            const float w2=omega*omega;
            sinc=1.0f+w2*(-1.0f/6.0f+w2/120.0f);
        } else {
            sinc=sin(phase)/omega;
        }
        total-=sinc;
    }
    return total;
}

inline Result evaluate(Point p,constant float2* nodes,uint order) {
    const float pi=3.14159265358979323846f;
    const float nf=float(p.n);
    const float parity=p.n&1u?-1.0f:1.0f;
    const float lambdaPower=powr(p.lambda,p.p);
    const float ratio=2.0f*p.lambda*lambdaPower;
    float length=1.0f-2.0f*p.lambda;
    if (p.ternary!=0u) length=1.0f/3.0f;
    float s=2.0f*pi*pi*nf*nf*p.d*powr(0.5f*length,p.p);
    float volume=length,product=1.0f,scale=1.0f;
    float coarse=0,fine=0,tail=INFINITY,approximation=0;
    uint den=3u,levels=0;
    for (uint j=0;j<96u;++j) {
        const float omega=pi*nf*length;
        const bool modular=p.ternary!=0u&&den!=0u;
        const float phase=modular?pi*float(p.n%(2u*den))/float(den):omega;
        const float fineError=s*s/(2.0f*(2.0f*p.p+1.0f))
                   +s*omega*omega/((p.p+1.0f)*(p.p+2.0f)*(p.p+3.0f));
        const float weight=volume*product*parity;
        float value;
        if (abs(weight)*fineError<=2e-14f) {
            value=-s/(p.p+1.0f);
            approximation+=abs(weight)*fineError;
        } else {
            value=gap_integral(nodes,order,phase,omega,s,p.p);
        }
        const bool geometricCoarse=p.ternary!=0u?modular&&p.n>=den:nf*length>=1.0f;
        if (geometricCoarse||s>=1.0f) coarse+=weight*value;
        else fine+=weight*value;
        levels=j+1u;
        tail=volume*s*ratio/((1.0f-ratio)*(p.p+1.0f));
        if (tail<=1e-13f) break;
        const float center=modular?2.0f*pi*float(p.n%den)/float(den):pi*nf*(1.0f-p.lambda)*scale;
        product*=cos(center);
        volume*=2.0f*p.lambda;
        length*=p.lambda;
        scale*=p.lambda;
        s*=lambdaPower;
        den=den!=0u&&den<1000000u?3u*den:0u;
    }
    return Result{coarse+fine,coarse,fine,tail,approximation,levels};
}

kernel void cantor(constant float2* nodes [[buffer(0)]],device Row* output [[buffer(1)]],
                   constant uint4& params [[buffer(2)]],
                   uint id [[thread_position_in_grid]],uint lane [[thread_index_in_simdgroup]]) {
    const uint i=params.z+id;
    const uint r=hash32(i^0xe746a3f5u^params.w);
    const uint powerSlot=i%11u;
    const float powers[11]={0.5f,0.75f,1,1.25f,1.5f,2,2.5f,3,4,6,8};
    const uint ternary=(i/11u)%3u==0u;
    const float lam=ternary!=0u?1.0f/3.0f:mix(0.23f,0.46f,unit(r,1));
    const float d=exp2(mix(-10.0f,5.0f,unit(r,2)));
    const uint pattern=(i/33u)%4u;
    const uint ns[6]={3u,9u,27u,81u,243u,729u};
    const uint pow3=ns[hash32(r+697921u)%6u];
    uint n=1u+hash32(r+10931u)%1024u;
    if (pattern==0u) n=pow3;
    if (pattern==1u) n=(hash32(r+216997u)&1u) != 0u ? pow3+1u : pow3-1u;
    const Point p=Point{lam,powers[powerSlot],d,n,ternary,pattern};
    const Result y=evaluate(p,nodes,params.y);
    const bool good=isfinite(y.full)&&isfinite(y.tail)&&isfinite(y.approximation);
    const uint group=id/32u;
    for (uint j=0;j<10u;++j) {
        const bool cancel=y.coarse*y.fine<0 && min(abs(y.coarse),abs(y.fine))>=1e-5f
                            && abs(y.full)<min(abs(y.coarse),abs(y.fine));
        const bool enabled=j==0u || (j==1u&&ternary!=0u&&pattern==0u&&n>=27u)
            || (j==2u&&ternary!=0u&&pattern==1u&&n>=26u)
            || (j==3u&&p.p<2.0f&&n>=128u) || (j==4u&&p.p==2.0f&&n>=128u)
            || (j==5u&&p.p>2.0f&&n>=128u) || (j==6u&&abs(y.fine)>abs(y.coarse))
            || (j==7u&&cancel) || (j==8u&&lam>=0.42f&&n>=128u) || j==9u;
        float metric=abs(y.full);
        if (j==6u) metric=abs(y.fine);
        if (j==7u) metric=min(abs(y.coarse),abs(y.fine));
        if (j==9u) metric=float(y.levels);
        const float val=enabled&&good?metric:-INFINITY;
        const float mx=simd_max(val);
        const uint winner=simd_min(val==mx?lane:32u);
        const uint count=simd_sum(uint(enabled));
        const uint bad=simd_sum(uint(enabled&&!good));
        if (lane==winner) output[group*10u+j]=Row{p,y.full,y.coarse,y.fine,y.tail,y.approximation,mx,y.levels,count,bad,i,params.w};
    }
}
