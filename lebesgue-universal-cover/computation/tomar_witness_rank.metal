#include <metal_stdlib>
using namespace metal;
constant uint directions=64;
// Reuleaux support material for n odd. For normals in a vertex's normal cone
// choose the vertex; between adjacent cones choose the unit arc about the
// vertex nearest to the opposite normal.
float2 contact(int n, float theta, float phi, float2 normal) {
    constexpr float tau=6.2831853071795864769f;
    float spacing=tau/float(n),radius=0.5f/cos(0.25f*spacing);
    float delta=phi-theta;
    float v=theta+spacing*rint(delta/spacing);
    if(abs(phi-v)<=0.25f*spacing) return radius*float2(cos(v),sin(v));
    v=theta+spacing*rint((delta+0.5f*tau)/spacing);
    return normal+radius*float2(cos(v),sin(v));
}
struct Site {float2 point,material; uint owner;};
kernel void rankWitness(device const float* poses [[buffer(0)]],
                        device const int* moves [[buffer(1)]],
                        device float* out [[buffer(2)]],
                        constant uint& nMoves [[buffer(3)]],
                        uint2 group [[threadgroup_position_in_grid]],
                        uint at [[thread_index_in_threadgroup]]) {
    const uint row=group.x, move=group.y;
    device const float* center=poses+22*row;
    device const float* width=center+11;
    device const int* direction=moves+11*move;
    constexpr float tau=6.2831853071795864769f;
    float phi=tau*(float(at)+0.31783797f)/float(directions);
    float2 u=float2(cos(phi),sin(phi));
    Site best={.5f*u,.5f*u,0};
    float h=.5f;
    for(uint b=1;b<=4;b++){
        int j=3*int(b-1),n=2*int(b)+1;
        float theta=(b==1)?0.f:center[j-1];
        float angle=(b==1)?0.f:width[j-1]*direction[j-1];
        float2 t=float2(center[j],center[j+1]), shift=float2(width[j]*direction[j],width[j+1]*direction[j+1]);
        float2 a=contact(n,theta+angle,phi,u);
        float val=dot(t+shift+a,u);
        if(val>h){
            h=val;
            float ca=cos(angle),sa=sin(angle);
            float2 ref=float2(ca*a.x+sa*a.y,-sa*a.x+ca*a.y);
            best={t+ref,ref,b};
        }
    }
    threadgroup Site chain[directions];
    chain[at]=best;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    // All arithmetic here ranks candidate material chains only.
    if(at!=0)return;
    float twice=0.f,grad[11]={};
    for(uint k=0;k<directions;k++){
        Site a=chain[k], next=chain[(k+1)%directions],prev=chain[(k+directions-1)%directions];
        twice+=a.point.x*next.point.y-a.point.y*next.point.x;
        if(a.owner==0) continue;
        int j=3*int(a.owner-1);
        float2 d=next.point-prev.point;
        grad[j]+=d.y;grad[j+1]-=d.x;
        if(a.owner>1)grad[j-1]-=dot(a.material,d);
    }
    for(uint j=0;j<11;j++)twice-=abs(grad[j])*width[j];
    out[row*nMoves+move]=0.5f*twice;
}
