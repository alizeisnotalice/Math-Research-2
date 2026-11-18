#include <metal_stdlib>
using namespace metal;

struct Hit {
    float a,b,p,q,score,aux;
    uint caseId,mode,index,tested,negative,nonfinite;
};
constant float norm=0.3989422804014327f;
inline uint hash32(uint x) {
    x^=x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float unit(uint x) {return float(hash32(x)>>8)*0x1p-24f;}
inline float em1(float x) {
    if (abs(x)>.1f) return exp(x)-1.0f;
    return x*(1.0f+x*(.5f+x*(1.0f/6.0f+x*(1.0f/24.0f+x*(1.0f/120.0f+x/720.0f)))));
}
inline float gaussian(float u,float v) {return norm*rsqrt(v)*exp(-0.5f*u*u/v);}
inline float keep_heldout(float t) {
    return abs(round(t*8.0f)-t*8.0f)<0.0001f ? t+0.0001f : t;
}

inline float2 transform(device const float* off,uint n,float R,float dx,float a,float b,
                        constant float2* nw,uint K) {
    float2 accum=0.0f;
    const float real=-.5f*(a*a-b*b),imag=-a*b;
    for(uint i=0;i+1u<n;++i) {
        const float lo=-R+dx*float(i),v0=off[i],dv=off[i+1]-v0;
        if (v0==0.0f&&dv==0.0f) continue;
        for(uint j=0;j<K;++j) {
            const float u=nw[j].x,x=lo+dx*u,v=v0+dv*u,d=x-a;
            const float w=norm*dx*nw[j].y*exp(-.5f*d*d),f=exp(real*v);
            if (b==0.0f) {
                accum.x+=w*em1(real*v);
            } else {
                const float h=sin(0.5f*imag*v),g=em1(real*v)-2.0f*f*h*h;
                const float s=f*sin(imag*v),phase=b*d;
                const float c=cos(phase),sz=sin(phase);
                accum+=w*float2(c*g-sz*s,sz*g+c*s);
            }
        }
    }
    return accum;
}

inline float Qsmooth(device const float* off,uint n,float R,float dx,float t,float y,float eps,
                     constant float2* nw,uint K) {
    float acc=t<1.0f?gaussian(y,t+eps):0.0f;
    for(uint i=0;i+1u<n;++i) {
        const float lo=-R+dx*float(i),v0=1.0f+off[i],v1=1.0f+off[i+1],dv=v1-v0;
        if (v0==1.0f&&v1==1.0f) continue;
        float l=0.0f,h=1.0f;
        if (dv>0.0f) h=clamp((t-v0)/dv,0.0f,1.0f);
        if (dv<0.0f) l=clamp((t-v0)/dv,0.0f,1.0f);
        if (dv==0.0f&&v0>t) h=0.0f;
        if (t>=1.0f) {
            for(uint j=0;j<K;++j) {
                const float x=lo+dx*nw[j].x;
                acc+=norm*dx*nw[j].y*exp(-.5f*x*x)*gaussian(y-x,t-1.0f+eps);
            }
        }
        if (h<=l) continue;
        for(uint j=0;j<K;++j) {
            const float u=l+(h-l)*nw[j].x,x=lo+dx*u,v=v0+dv*u;
            acc-=norm*dx*(h-l)*nw[j].y*exp(-.5f*x*x)*gaussian(y-x,max(t-v,0.0f)+eps);
        }
    }
    return acc;
}

inline void append(thread float2& acc,float value) {
    if (value>acc.x) {
        acc.y=1.0f+acc.y*exp(acc.x-value);
        acc.x=value;
    } else if (isfinite(value)) acc.y+=exp(value-acc.x);
}
inline float logem(float v) {
    if (v>=18.0f) return 0.0f;
    if (v<=0.0f) return -INFINITY;
    return log(-em1(-v));
}
inline float2 tailLogs(device const float* off,uint n,float R,float dx,float t,
                       constant float2* nw,uint K) {
    float2 plus=float2(-INFINITY,0),minus=plus;
    const float q=t*t/2;
    for (uint i=0;i+1u<n;++i) {
        const float l=-R+dx*float(i),r=l+dx,b0=off[i],b1=off[i+1],sl=(b1-b0)/dx;
        if (b0==0.0f&&b1==0.0f) continue;
        for (uint m=0;m<2u;++m) {
            const float side=m==0u?-1.0f:1.0f;
            if (max(side*b0,side*b1)<=0.0f) continue;
            float a=l,b=r;
            if (side*b0<0.0f) a=l-b0/sl;
            if (side*b1<0.0f) b=l-b0/sl;
            const float effective=t-(m==0u?q*sl:0.0f);
            const float xpeak=clamp(effective,a,b);
            const float vmax=t*xpeak-0.5f*xpeak*xpeak+(m==0u?-q*(b0+sl*(xpeak-l)):0.0f);
            for(uint k=0;k<2u;++k) {
                const float dir=k==0u?-1.0f:1.0f;
                float width=k==0u?xpeak-a:b-xpeak;
                if (width<=0.0f) continue;
                const float speed=max(0.0f,-dir*(effective-xpeak));
                if (speed*width>42.0f) width=42.0f/speed;
                const float factor=vmax+log(norm*width);
                for(uint j=0;j<K;++j) {
                    const float u=dir*width*nw[j].x,x=xpeak+u;
                    const float w=b0+sl*(x-l);
                    const float val=factor+(effective-xpeak)*u-0.5f*u*u
                                    +logem(q*side*w)+log(nw[j].y);
                    if (m==0u) append(plus,val);
                    else append(minus,val);
                }
            }
        }
    }
    return float2(plus.x+log(plus.y),minus.x+log(minus.y));
}

kernel void validation(device const float* profiles [[buffer(0)]],constant float2* nw [[buffer(1)]],
                  device Hit* out [[buffer(2)]],constant uint4& P [[buffer(3)]],
                  uint id [[thread_position_in_grid]],uint lane [[thread_index_in_simdgroup]]) {
    const uint groups=9u*P.x, runs=P.w;
    const uint index=P.y+id,gid=id/32u;
    const uint stratum=(gid/runs)%groups, which=stratum/9u,mode=stratum%9u;
    device const float* row=profiles+which*P.z;
    const uint knots=uint(row[2]),K=8u;
    const float R=row[0],dx=row[1];
    device const float* off=row+5;
    const uint seed=hash32(index+8349551u);
    const float u=unit(seed),v=unit(seed+12372801u),w=unit(seed+43886493u);
    float a=0.0f,b=0.0f,p=0.0f,q=0.0f,score=0.0f,aux=0.0f;
    bool neg=false;
    if (mode<=3u) {
        if (mode<=1u) p=keep_heldout(mix(.10f,8.0f,u))*(mode==0u?1.0f:-1.0f);
        if (mode==2u) q=mix(.10f,8.0f,u)*(v<.5f?1.0f:-1.0f);
        if (mode==3u) {
            p=mix(.15f,4.0f,u)*(w<.5f?1.0f:-1.0f);
            q=mix(.10f,4.0f,v);
        }
        const float2 res=transform(off,knots,R,dx,p,q,nw,K);
        a=res.x; b=res.y; score=length(res); neg=a<0.0f;
    } else if (mode==6u) {
        p=mix(.52f,1.53f,u); q=mix(-4.0f,4.0f,v);
        aux=w<.5f?.03125f:.125f;
        a=Qsmooth(off,knots,R,dx,p,q,aux,nw,K);
        score=-a; neg=a<0.0f;
    } else {
        const float t=mode<=5u?mix(8.0f,32.0f,u):16.0f*(R+.125f)+mix(16.0f,64.0f,u);
        p=t*((mode==4u||mode==7u)?1.0f:-1.0f);
        const float2 logs=tailLogs(off,knots,R,dx,p,nw,K);
        // Terms can be absent in a sign-definite graph, retain a finite sentinel.
        a=max(-1e20f,logs.x); b=max(-1e20f,logs.y);
        neg=a<b;
        score=mode<=5u?max(a,b)-0.5f*p*p:a-b;
        aux=a-b;
    }
    const bool good=isfinite(a)&&isfinite(b)&&isfinite(score);
    const float value=good?score:-INFINITY;
    const float mx=simd_max(value);
    const uint winner=simd_min(value==mx?lane:32u);
    const uint bad=simd_sum(uint(!good)),negCount=simd_sum(uint(neg&&good));
    if (lane==winner) out[gid]=Hit{a,b,p,q,mx,aux,which,mode,index,32u,negCount,bad};
}
