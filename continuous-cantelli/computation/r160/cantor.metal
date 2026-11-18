#include <metal_stdlib>
using namespace metal;

inline float3 moment(float y,float x,float iv,float normalized) {
    const float z=x-y,k=normalized*exp(-0.5f*z*z*iv),v=z*iv;
    return float3(k,-v*k,(v*v-iv)*k);
}
kernel void cantor(device const float4* queries [[buffer(0)]],
                   device const float2* gaps [[buffer(1)]],
                   constant float2* nodes [[buffer(2)]],
                   device float* result [[buffer(3)]],
                   uint id [[thread_position_in_grid]]) {
    const float4 row=queries[id];
    const float x=row.x,iv=1.0f/row.y,norm=.3989422804014327f*sqrt(iv);
    const uint stage=uint(row.z),order=uint(row.w),count=(1u<<stage)-1u;
    const uint boundary=count==0u?0u:(count-1u)/2u;
    float3 all=0.0f,gap=0.0f,current=0.0f,error=0.0f;
    for(uint j=0;j<16u;++j)
        all+=nodes[j].y*moment(nodes[j].x,x,iv,norm);
    for(uint i=0;i<count;++i) {
        const float2 g=gaps[i];
        float3 sum=0.0f;
        if(order==1u) {
            sum=g.y*moment(g.x,x,iv,norm);
        } else {
            for(uint j=16u;j<20u;++j)
                sum+=g.y*nodes[j].y*moment(g.x+g.y*nodes[j].x,x,iv,norm);
        }
        if(i<boundary) {
            const float3 temp=sum-error,step=gap+temp;
            error=(step-gap)-temp;
            gap=step;
        } else {
            current+=sum;
        }
    }
    const float3 surface=all-gap-current;
    result[id*6u]=surface.x;
    result[id*6u+1u]=surface.y;
    result[id*6u+2u]=surface.z;
    result[id*6u+3u]=gap.x+current.x;
    result[id*6u+4u]=gap.z+current.z;
    result[id*6u+5u]=current.x;
}
