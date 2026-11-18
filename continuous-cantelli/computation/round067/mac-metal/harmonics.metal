#include <metal_stdlib>
using namespace metal;
struct TablePoint { float4 c1,c2,s1,s2; };
struct Settings { uint panels,family,seed; float amplitude; };
struct Candidate { float4 metrics,c1,c2,s1,s2; };

uint mixBits(uint x) {
    x ^= x>>16;
    x *= 0x7feb352du;
    x ^= x>>15;
    x *= 0x846ca68bu;
    return x^(x>>16);
}
float signedDraw(thread uint& state) {
    state=mixBits(state+0x9e3779b9u);
    return (float(state>>8)/8388608.0f)-1.0f;
}

kernel void scan_harmonics(device const TablePoint* table [[buffer(0)]],
                           constant Settings& settings [[buffer(1)]],
                           device Candidate* output [[buffer(2)]],
                           uint idx [[thread_position_in_grid]]) {
    uint state=mixBits(idx+104729u*settings.seed);
    float c[8]={0},s[8]={0};
    uint f=settings.family;
    float norm2=0.0f;
    for(uint k=0;k<8;++k) {
        bool active=(f==0 && k<3) || ((f==1 || f==2 || f==4) && (k==3 || k==4 || k==6)) || (f==3 && (k==3 || k==7));
        if(active) {
            c[k]=signedDraw(state);
            s[k]=signedDraw(state);
            norm2+=c[k]*c[k]+s[k]*s[k];
        }
    }
    float scale=settings.amplitude/sqrt(max(norm2,1.e-15f));
    for(uint k=0;k<8;++k) {c[k]*=scale;s[k]*=scale;}
    if(f==2) {
        for(uint n=1;n<=3;++n) {
            float re=0,im=0;
            for(uint k=3;k+n<8;++k) {
                uint j=k+n;
                re+=c[j]*c[k]+s[j]*s[k];
                im+=s[j]*c[k]-c[j]*s[k];
            }
            float halfRange=settings.amplitude*settings.amplitude*settings.amplitude*0.75f;
            c[n-1]=float(n*n)*re*0.25f+halfRange*signedDraw(state);
            s[n-1]=float(n*n)*im*0.25f+halfRange*signedDraw(state);
        }
    }
    if(f==4) {
        for(uint n=0;n<3;++n) {
            c[n]=0.3f*settings.amplitude*signedDraw(state);
            s[n]=0.3f*settings.amplitude*signedDraw(state);
        }
    }
    float4 c1=float4(c[0],c[1],c[2],c[3]),c2=float4(c[4],c[5],c[6],c[7]);
    float4 s1=float4(s[0],s[1],s[2],s[3]),s2=float4(s[4],s[5],s[6],s[7]);
    float3 x=0,y=0,z=0;
    for(uint j=0;j<settings.panels;++j) {
        TablePoint t=table[j];
        float a=dot(c1,t.c1)+dot(c2,t.c2)+dot(s1,t.s1)+dot(s2,t.s2);
        float3 mass=exp(-a*float3(0.5f,2.0f,4.5f));
        x+=mass*t.c1.xyz;
        y+=mass*t.s1.xyz;
        z+=mass;
    }
    float3 res=sqrt(x*x+y*y)/z;
    output[idx]={float4(res,max(max(res.x,res.y),res.z)),c1,c2,s1,s2};
}
