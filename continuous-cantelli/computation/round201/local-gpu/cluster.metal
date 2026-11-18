#include <metal_stdlib>
using namespace metal;

constant uint order=64u;
inline float start(uint j,uint m,float rho) {
    return j==0u?0.0f:.9f+float(j-1u)*(.1f-rho)/float(m-2u);
}
inline float J(float z,float k,uint ki,constant float2* nw) {
    const float zz=z*z;
    float acc=0.0f;
    for(uint j=0;j<order;++j) {
        const float2 node=nw[ki*order+j];
        const float x=zz*node.x;
        acc+=node.y*(1.0f-x)*exp(-.5f*x);
    }
    return acc/(2.0f*k);
}
inline float shell(float y,float D,float m,float logm,float k,float gamma,uint r) {
    const float a=k+.5f, t=2.0f/(D*D*y*y);
    const float small=t*exp(-4.0f*float(r)*logm);
    const float c1=(k+1.0f)*a,c2=.5f*(k+1.0f)*(k+2.0f)*a*(a+1.0f);
    return -gamma*powr(t,a)*exp(-(1.0f+(4.0f*a-1.0f)*float(r))*logm)
           *(1.0f+small*(c1+small*c2))+0.0f*m;
}
kernel void cluster(device const float4* config [[buffer(0)]],
                    constant float2* nw [[buffer(1)]],
                    device float* result [[buffer(2)]],
                    uint id [[thread_position_in_grid]]) {
    const float4 p=config[id];
    const uint M=uint(p.x),ki=uint(p.z),depth=uint(p.w);
    const float m=float(M),D=p.y,k=ki==0u?.25f:.5f,rho=1.0f/(m*m);
    const float mean=(m-1.0f)/m*(1.9f-rho)/(2.0f*(1.0f-rho));
    const float gamma=ki==0u?1.2254167024651776f:1.0f;
    float sum=0.0f;
    for(uint i=0;i<M;++i) {
        const float a=start(i,M,rho);
        if(depth==1u) {
            sum+=J(D*(a+rho*mean),k,ki,nw);
        } else {
            for(uint j=0;j<M;++j) {
                const float x=a+rho*(start(j,M,rho)+rho*mean);
                sum+=J(D*x,k,ki,nw);
            }
        }
    }
    sum/=depth==1u?m:m*m;
    result[6u*id]=sum;
    const float logm=log(m);
    for(uint r=1u;r<=8u;++r) {
        for(uint i=1u;i<M;++i) {
            const float a=start(i,M,rho);
            sum+=shell(a+rho*mean,D,m,logm,k,gamma,r);
        }
        if(r==1u) result[6u*id+1u]=sum;
        if(r==2u) result[6u*id+2u]=sum;
        if(r==4u) result[6u*id+3u]=sum;
        if(r==8u) result[6u*id+4u]=sum;
    }
    result[6u*id+5u]=J(0.0f,k,ki,nw);
}
