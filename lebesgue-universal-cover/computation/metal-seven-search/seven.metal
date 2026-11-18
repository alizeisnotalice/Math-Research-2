#include <metal_stdlib>
using namespace metal;
constant float pi=3.14159265358979323846f;
float reuleaux(float a, float n) {
    float period=2*pi/n;
    float z=abs(a-rint(a/period)*period);
    float r=0.5f/cos(pi/(2*n));
    return z<=pi/(2*n)?r*cos(z):1-r*cos(pi/n-z);
}
kernel void seven(device const float *poses [[buffer(0)]],constant float2 *cloud [[buffer(1)]],
                  device float *areas [[buffer(2)]],uint t [[thread_index_in_threadgroup]],
                  uint3 group [[threadgroup_position_in_grid]]) {
    device const float *p=poses+18*group.x;
    threadgroup float h[1024], sums[128];
    float2 rp=float2(cos(p[11]),sin(p[11])),rq=float2(cos(p[14]),sin(p[14]));
    for(uint j=t;j<1024;j+=128) {
        float a=2*pi*j/1024; float2 u=float2(cos(a),sin(a));
        float mx=max(0.5f,reuleaux(a,3)+dot(u,float2(p[0],p[1])));
        for(int k=0;k<3;k++) {
            int i=2+3*k;
            mx=max(mx,reuleaux(a-p[i],5+2*k)+dot(u,float2(p[i+1],p[i+2])));
        }
        float2 up=float2(dot(u,rp),u.y*rp.x-u.x*rp.y);
        float2 uq=float2(dot(u,rq),(u.y*rq.x-u.x*rq.y)*p[17]);
        float hp=-2,hq=-2;
        for(int i=0;i<30;i++) hp=max(hp,dot(up,cloud[i]));
        for(int i=30;i<41;i++) hq=max(hq,dot(uq,cloud[i]));
        h[j]=max(mx,max(hp+dot(u,float2(p[12],p[13])),hq+dot(u,float2(p[15],p[16]))));
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float localArea=0;
    for(uint j=t;j<1024;j+=128) {
        float v=h[j],diff=v-h[(j+1)%1024];
        localArea+=v*v*tan(pi/1024)-diff*diff/(2*sin(2*pi/1024));
    }
    sums[t]=localArea;
    for(uint step=64;step;step/=2) {
        threadgroup_barrier(mem_flags::mem_threadgroup);
        if(t<step)sums[t]+=sums[t+step];
    }
    if(t==0)areas[group.x]=sums[0];
}
uint shuffle_bits(uint x) {
    x^=x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
float rnd(uint x){return float(shuffle_bits(x)>>8)*(1.0f/16777216);}
kernel void mutate(device const float* centers [[buffer(0)]],device float* poses [[buffer(1)]],
                   constant uint2 &params [[buffer(2)]],uint id [[thread_position_in_grid]]) {
    uint child=id/18,k=id%18; float value=centers[(child/512)*18+k];
    uint base=shuffle_bits(child+params.x*8192u)^params.y;
    if(k<17 && child%512 && rnd(base+37u*k)<.4f){
        float sigma=.008f*pow(.96f,float(params.x%96))*pow(10.f,1.5f*(rnd(base+999)-.5f));
        float z=sqrt(-2*log(max(1.e-7f,rnd(base+91*k+11))))*cos(6.2831853f*rnd(base+53*k+13));
        value+=z*sigma*((k==2||k==5||k==8||k==11||k==14)?15.f:1.f);
    }
    poses[id]=value;
}
kernel void select_best(device const float* poses [[buffer(0)]],device const float* area [[buffer(1)]],
                        device float* centers [[buffer(2)]],device float* scores [[buffer(3)]],
                        uint j [[thread_index_in_threadgroup]],uint3 g [[threadgroup_position_in_grid]]){
    threadgroup uint winner[128];uint row=g.x*512+j;
    for(uint k=1;k<4;k++){uint next=g.x*512+j+128*k;if(area[next]<area[row])row=next;}
    winner[j]=row;
    for(uint step=64;step;step/=2){threadgroup_barrier(mem_flags::mem_threadgroup);
       if(j<step && area[winner[j+step]]<area[winner[j]])winner[j]=winner[j+step];}
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if(j<18)centers[g.x*18+j]=poses[winner[0]*18+j];if(j==0)scores[g.x]=area[winner[0]];
}
