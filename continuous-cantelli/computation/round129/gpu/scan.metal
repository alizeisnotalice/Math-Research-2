#include <metal_stdlib>
using namespace metal;
struct Hit {
    float a,b,p,q,score,aux,ap,am,pp,pm;
    uint caseId,mode,index,tested,negative,nonfinite;
};
constant float norm=.3989422804014327f;
inline uint hash32(uint x) {
    x^=x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float unit(uint x) {return float(hash32(x)>>8)*0x1p-24f;}
inline float em1(float x) {
    if (abs(x)>.1f) return exp(x)-1.0f;
    return x*(1.0f+x*(.5f+x*(1.0f/6.0f+x*(1.0f/24.0f+x*(1.0f/120.0f+x/720.0f)))));
}
inline float offset(float r,float R,float a,float p) {
    const float u=R/r;
    return a*powr(u,p)*(1.0f-u);
}
inline void residual(thread float2& accum,float x,float v,float weight,float a,float b) {
    const float d=x-a,w=norm*weight*exp(-.5f*d*d),rho=-.5f*(a*a-b*b)*v;
    if (b==0.0f) {
        accum.x+=w*em1(rho);
    } else {
        const float theta=-a*b*v,phase=b*d;
        const float sn=sin(.5f*theta),e=exp(rho);
        const float g=em1(rho)-2.0f*e*sn*sn,s=e*sin(theta);
        const float c=cos(phase),s2=sin(phase);
        accum+=w*float2(c*g-s2*s,s2*g+c*s);
    }
}
inline float2 compact(device const float* off,uint n,float R,float dx,float a,float b,
                     constant float2* nw) {
    float2 res=0.0f;
    for(uint i=0;i+1u<n;++i) {
        const float x0=-R+dx*float(i),v=off[i],dv=off[i+1]-v;
        for(uint j=0;j<8u;++j)
            residual(res,x0+dx*nw[j].x,v+dv*nw[j].x,dx*nw[j].y,a,b);
    }
    return res;
}
inline void append(thread float2& acc,float value) {
    if (value>acc.x) {
        acc.y=1.0f+acc.y*exp(acc.x-value);
        acc.x=value;
    } else if (isfinite(value)) acc.y+=exp(value-acc.x);
}
inline float whole_log(device const float* off,uint n,float R,float dx,float t,
                       float a,float p,constant float2* nw) {
    float2 sum=float2(-INFINITY,0.0f);
    const float q=.5f*t*t,mag=abs(t),dir=t>0.0f?1.0f:-1.0f;
    if(mag<64.0f) for(uint i=0;i+1u<n;++i) {
        const float x0=-R+dx*float(i),v=off[i],dv=off[i+1]-v;
        for(uint j=0;j<8u;++j) {
            const float x=x0+dx*nw[j].x,z=x-t;
            append(sum,log(norm*dx*nw[j].y)-.5f*z*z-q*(v+dv*nw[j].x));
        }
    }
    const float l=max(R-mag,-20.0f),w=(20.0f-l)/20.0f;
    for(uint i=0;i<20u;++i) {
        for(uint j=0;j<8u;++j) {
            const float z=l+w*(float(i)+nw[j].x);
            const float b=offset(mag+z,R,a,p);
            append(sum,log(norm*w*nw[j].y)-.5f*z*z-q*b);
        }
    }
    return sum.x+log(sum.y)+0.0f*dir;
}
kernel void validation(device const float* profiles [[buffer(0)]],constant float2* nw [[buffer(1)]],
                  device Hit* out [[buffer(2)]],constant uint4& P [[buffer(3)]],
                  uint id [[thread_position_in_grid]],uint lane [[thread_index_in_simdgroup]]) {
    const uint index=P.y+id,gid=id/32u;
    const uint stratum=(gid/P.w)%(9u*P.x),which=stratum/9u,mode=stratum%9u;
    device const float* row=profiles+which*P.z;
    const uint knots=uint(row[2]);
    const float R=row[0],dx=row[1];
    const uint seed=hash32(index+213952837u);
    const float u=unit(seed),v=unit(seed+6280183u),w=unit(seed+4188973u);
    const float ap=(v<.5f?1.0f:-1.0f)*mix(.025f,.22f,unit(seed+4107231u));
    const float am=(w<.5f?1.0f:-1.0f)*mix(.025f,.22f,unit(seed+7542609u));
    const float pl=(mode==4u||mode==5u)?.55f:(mode==8u?2.1f:.55f);
    const float ph=(mode==4u||mode==5u)?1.9f:3.5f;
    const float pp=(mode==6u||mode==7u)?2.0f:mix(pl,ph,unit(seed+1598269u));
    const float pm=(mode==6u||mode==7u)?2.0f:mix(pl,ph,unit(seed+2444819u));
    float p=0.0f,q=0.0f,a=0.0f,b=0.0f,aux=0.0f,score=0.0f;
    if(mode<4u) {
        if(mode<2u) {
            p=mix(.11f,7.98f,u);
            if(abs(round(p*8.0f)-p*8.0f)<.0001f) p+=.0001f;
            if(mode==1u) p=-p;
        }
        if(mode==2u) q=mix(.19f,7.95f,u);
        if(mode==3u) {
            p=mix(.17f,3.95f,u)*(w<.5f?1.0f:-1.0f);
            q=mix(.12f,3.95f,v);
        }
        float2 z=compact(row+5,knots,R,dx,p,q,nw);
        for(uint k=0;k<2u;++k) for(uint i=0;i<12u;++i) for(uint j=0;j<8u;++j) {
            const float r=R+float(i)+nw[j].x,sg=k==0u?1.0f:-1.0f;
            residual(z,sg*r,offset(r,R,k==0u?ap:am,k==0u?pp:pm),nw[j].y,p,q);
        }
        a=z.x;b=z.y;score=length(z);
    } else {
        const bool pos=mode==4u||mode==6u||(mode==8u&&v<.5f);
        p=mix(22.0f,112.0f,u)*(pos?1.0f:-1.0f);
        const float A=pos?ap:am,B=pos?pp:pm;
        a=whole_log(row+5,knots,R,dx,p,A,B,nw);
        b=-.5f*p*p*offset(abs(p),R,A,B);
        aux=-.5f*A*powr(R,B);
        score=abs(a);
    }
    const bool good=isfinite(a)&&isfinite(b)&&isfinite(score);
    const float value=good?score:-INFINITY;
    const float mx=simd_max(value);
    const uint winner=simd_min(value==mx?lane:32u);
    const uint bad=simd_sum(uint(!good)),negs=simd_sum(uint(a<0.0f&&good));
    if (lane==winner) out[gid]=Hit{a,b,p,q,mx,aux,ap,am,pp,pm,which,mode,index,32u,negs,bad};
}
