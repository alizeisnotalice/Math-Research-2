#include <metal_stdlib>
using namespace metal;
struct Core {float b,t,beta,Ccrit;};
struct Settings {uint panels,lambdas,offsets;};
float em1(float x){
 if(fabs(x)>0.08f)return exp(x)-1.;
 return x*(1.+x*(.5+x*(.1666666667+x*(.04166666667+x*(.00833333333+x*.00138888889)))));
}
float lp1(float x){
 if(fabs(x)>0.05f)return log(1.+x);
 return x*(1.+x*(-.5+x*(.3333333333+x*(-.25+x*(.2-x/6.)))));
}
float sinhcm1(float x){
 if(fabs(x)>.2)return sinh(x)/x-1.;
 float q=x*x;return q*(1./6.+q*(1./120.+q*(1./5040.+q/362880.)));
}
kernel void scan_laplace(device const Core *grid [[buffer(0)]],
                        constant Settings& cfg [[buffer(1)]],
                        device float4 *out [[buffer(2)]],
                        uint idx [[thread_position_in_grid]]) {
 uint point=cfg.lambdas*cfg.offsets;
 Core p=grid[idx/point];
 uint m=idx%point;
 float L=mix(log(32.),log(65536.),float(m/cfg.offsets)/float(cfg.lambdas-1));
 float lambda=exp(L);
 float off=2.*float(m%cfg.offsets)/float(cfg.offsets-1)-1.;
 float c=p.Ccrit+off;
 float K=p.beta*L+c;
 float tail=0,normal=0,cap=0;
 for(uint k=0;k<cfg.panels;k++){
    float f=(float(k)+.5)/float(cfg.panels);
    float y=16.*f-8.;
    float z=y-p.b,q=z/lambda;
    float inv=1./(1.+q);
    float residual=p.b*z*z*inv/lambda-.5*(p.beta*lp1(q)-K*q*(2.+q))*inv*inv;
    float weight=exp(-.5*y*y);
    normal+=weight;
    tail+=weight*em1(residual);
    float s=64.*f, v=s/lambda;
    float e=-.5*v*v/p.t, a=p.b*v/p.t;
    cap+=s*exp(-s)*(em1(e)+exp(e)*sinhcm1(a));
 }
 float logTail=lp1(tail/normal);
 float logCap=lp1(64.*cap/float(cfg.panels));
 float margin=(2.-.5*p.beta)*L-.5*off+logTail-logCap;
 float first=p.b*(p.b*p.b+1.+.5*p.beta-K);
 out[idx]=float4(margin,logTail,logCap,lambda*logTail-first);
}
