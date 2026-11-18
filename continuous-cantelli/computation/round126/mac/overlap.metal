#include <metal_stdlib>
using namespace metal;

struct Settings {uint seed, panels, mode, count;};
struct Candidate {float4 geometry,shape,overlap,values;};
uint mixbits(uint x) {
    x^=x>>16;x*=0x7feb352du;x^=x>>15;x*=0x846ca68bu;return x^(x>>16);
}
float unit(uint x) {return float(mixbits(x)>>8)*0x1p-24f;}
float expm1Small(float x) {
    if (fabs(x)>.045f) return exp(x)-1;
    return x*(1+x*(.5f+x*(1.f/6+x*(1.f/24+x/120))));
}
float log1pSmall(float x) {
    if (fabs(x)>.04f) return log(1+x);
    return x*(1+x*(-.5f+x*(1.f/3+x*(-.25f+x/5))));
}

kernel void overlaps(constant Settings& cfg [[buffer(0)]],
                     device Candidate *out [[buffer(1)]],
                     uint idx [[thread_position_in_grid]]) {
    uint seed=mixbits(idx+cfg.seed)+0x9e3779b9u;
    float A=exp2(4.f+10.f*unit(seed));
    if (cfg.mode==1) A=exp2(5.f+5.f*unit(seed));
    if (cfg.mode==2) A=exp2(8.f+6.f*unit(seed));
    float mid=2.f*(unit(seed+1)-.5f);
    float sep=.12f+3.3f*unit(seed+2);
    float widthL=(.15f+.75f*unit(seed+3))*sep;
    float widthR=(.15f+.75f*unit(seed+4))*sep;
    float sigma=exp2(-1.5f+3.5f*unit(seed+5));
    float depth=4.f*(unit(seed+6)-.5f);
    float p=.15f+.82f*unit(seed+7);
    Candidate record={float4(A,mid,sep,widthL),float4(widthR,sigma,depth,p),float4(0),float4(-1)};
    out[idx]=record;
    float sqrtA=sqrt(A), invT=exp(-.5f*A);
    float4 I=float4(0);
    for (uint j=0;j<2;j++) {
       float center=mid+(j==0?-sep:sep);
       float width=j==0?widthL:widthR;
       float off=j==0?0.f:depth;
       for (uint k=0;k<cfg.panels;k++) {
          float y=width*(2.f*(float(k)+.5f)/float(cfg.panels)-1.f);
          float z=center+y;
          float x=(A+sqrtA*z)*invT;
          float rho=(y/sigma)*(y/sigma)+off;
          float weight=exp(-.5f*(x*x+rho))*width;
          float v=(2.f*sqrtA)*z+rho;
          float d=-z-rho/sqrtA;
          I += weight*float4(1.f,v,d,d*d);
       }
    }
    float norm=I.x, c=log(A)+2.f*log(.7978845608f*norm/(p*float(cfg.panels)));
    if (!isfinite(c)) return;
    for (uint j=0;j<2;j++) {
        float center=mid+(j==0?-sep:sep);
        float width=j==0?widthL:widthR;
        float y=clamp(-sqrtA*sigma*sigma,-width,width);
        float rho=(y/sigma)*(y/sigma)+(j==0?0.f:depth);
        if (A+2.f*sqrtA*(center+y)+rho+c<=0.f) return;
    }
    float meanD=I.z/norm-c/sqrtA;
    float meanD2=I.w/norm-2.f*c/sqrtA*(I.z/norm)+c*c/A;
    float meanF=1.f+(I.y/norm+c)/A;
    float variance=p*meanD2+p*p/(1.f-p)*meanD*meanD;
    float energy=p*meanF;
    const float steps[8]={.125f,.25f,.40f,.60f,.85f,1.15f,1.60f,2.25f};
    float bestRelative=-1.e30f, bestEnergy=-1.e30f;
    float bestLambda=0, bestQpos=0, bestQneg=0, finiteSigns=0;
    for (uint si=0;si<8;si++) {
       float a=steps[si], halfSquare=.5f*a*a;
       float2 m=float2(0);
       for (uint j=0;j<2;j++) {
          float center=mid+(j==0?-sep:sep);
          float width=j==0?widthL:widthR;
          float off=j==0?0.f:depth;
          for (uint k=0;k<cfg.panels;k++) {
             float y=width*(2.f*(float(k)+.5f)/float(cfg.panels)-1.f);
             float z=center+y;
             float x=(A+sqrtA*z)*invT;
             float rho=(y/sigma)*(y/sigma)+off;
             float weight=exp(-.5f*(x*x+rho))*width;
             float D=-z-(c+rho)/sqrtA;
             float F=1.f+(2.f*sqrtA*z+c+rho)/A;
             m+=weight*float2(expm1Small(a*D-halfSquare*F),expm1Small(-a*D-halfSquare*F));
          }
       }
       float2 ratio=(p/norm)*m;
       float2 mass=p+ratio;
       if (mass.x>.98f || mass.y>.98f) continue;
       const float rem=1.f-p;
       float gain=-rem*(log1pSmall(-ratio.x/rem)+log1pSmall(-ratio.y/rem))/(a*a);
       float lower=energy+gain;
       float benefit=lower-max(energy,variance);
       if (benefit>bestRelative) {
          bestRelative=benefit;bestEnergy=lower;bestLambda=a;bestQpos=mass.x;bestQneg=mass.y;
       }
       if (lower>1.f && energy<1.f && variance<1.f) finiteSigns+=1;
    }
    record.overlap=float4(bestLambda,bestQpos,bestQneg,bestEnergy);
    record.values=float4(energy,variance,bestRelative,finiteSigns);
    out[idx]=record;
}
