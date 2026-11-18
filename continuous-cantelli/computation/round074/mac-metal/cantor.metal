#include <metal_stdlib>
using namespace metal;
struct Settings { uint panels,base,seed,density; float power; };
struct Candidate {float4 scores;uint4 config;};
uint hash(uint x){x^=x>>16;x*=0x7feb352du;x^=x>>15;x*=0x846ca68bu;return x^(x>>16);}
float draw(thread uint& state){state=hash(state+0x9e3779b9u);return float(state>>8)/16777216.0f;}
float2 cmul(float2 a,float2 b){return float2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);}
float loss(float x){
    if(x>=0.05f)return 1.0f-exp(-x);
    return x*(1.0f+x*(-0.5f+x*(0.16666666667f+x*(-0.04166666667f+x*0.00833333333f))));
}
kernel void cantor_scan(constant Settings& cfg [[buffer(0)]],
                       device Candidate* result [[buffer(1)]],
                       uint idx [[thread_position_in_grid]]) {
    uint state=hash(idx+104729u*cfg.seed);
    uint B=cfg.base;
    float density=0.25f*(1.0f+float(cfg.density));
    uint mask=1u|(1u<<(B-1));
    for(uint d=1;d+1<B;++d) if(draw(state)<density)mask|=1u<<d;
    if(mask==((1u<<B)-1u))mask^=1u<<(B/2);
    uint depth=1u+min(9u,uint(draw(state)*10.0f));
    float height=exp(log(0.03f)+draw(state)*log(1000.0f));
    float ib=1.0f/float(B);
    uint count=popcount(mask);
    float largest=0.0f,first3=0.0f,rawmax=0.0f,firstRaw=0.0f;
    for(uint n=1;n<=6;++n){
        float2 val=0,prefix=float2(1.0f,0.0f);
        float mass=0.0f,angle=2.0f*M_PI_F*float(n),scale=1.0f,live=1.0f;
        float strength=0.5f*float(n*n)*height;
        for(uint k=0;k<depth;++k){
            uint left=0;
            for(uint right=1;right<B;++right){
                if(!(mask&(1u<<right)))continue;
                if(right>left+1){
                    float gap=float(right-left-1)*ib*scale;
                    float center=float(right+left+1)*0.5f*ib*scale;
                    float u=0.0f,v=0.0f;
                    float gain=strength*pow(gap,cfg.power);
                    for(uint j=0;j<cfg.panels;++j){
                        float z=(float(j)+0.5f)/(2.0f*float(cfg.panels));
                        float f=loss(gain*pow(z,cfg.power));
                        u+=f;
                        v+=f*cos(angle*gap*(z-0.5f));
                    }
                    float w=live*gap/float(cfg.panels)/scale;
                    mass+=w*u;
                    val+=w*v*cmul(prefix,float2(cos(angle*center),sin(angle*center)));
                }
                left=right;
            }
            if(k+1<depth){
                float2 m=0;
                for(uint d=0;d<B;++d)if(mask&(1u<<d)){
                    float t=angle*float(d)*ib*scale;
                    m+=float2(cos(t),sin(t));
                }
                prefix=cmul(prefix,m/float(count));
                live*=float(count)*ib;
                scale*=ib;
            }
        }
        float raw=length(val),score=raw/mass;
        if(n==1)firstRaw=raw;
        if(n<=3)first3=max(first3,score);
        rawmax=max(rawmax,raw);
        largest=max(largest,score);
    }
    result[idx]={float4(largest,first3,rawmax,firstRaw),uint4(mask,depth,as_type<uint>(height),0u)};
}
