#include <metal_stdlib>
using namespace metal;
struct Core {float b,t,share,Ccrit;};
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
 float c=p.Ccrit-2.*log(p.share);
 float K=4.*L+c;
 float logCapSlope=(p.b*p.b-p.Ccrit)*.5;
 float cEdge=exp(log(2.)-.5*log(2.*M_PI_F)-.5*p.b*p.b-log(1.-p.share)-logCapSlope);
 if (cEdge*lambda<16.) {
    out[idx]=float4(as_type<float>(0x7f7fffffu),0.,0.,0.);
    return;
 }
 float Dstar=-2.*(p.b*(p.b*p.b+3.-c)+(1.-p.share)/p.share*2./cEdge);
 float eta=8.*p.b+off;
 float tail=0,normal=0,cap=0,edge=0;
 for(uint k=0;k<cfg.panels;k++){
    float f=(float(k)+.5)/float(cfg.panels);
    float y=16.*f-8.;
    float z=y-p.b,q=z/lambda;
    float inv=1./(1.+q);
    float residual=p.b*z*z*inv/lambda-.5*(4.*lp1(q)-K*q*(2.+q))*inv*inv;
    residual+=.5*(eta*(L+lp1(q))+Dstar)*inv*inv*inv/lambda;
    float weight=exp(-.5*y*y);
    normal+=weight;
    tail+=weight*em1(residual);
    float s=64.*f, v=s/lambda;
    float e=-.5*v*v/p.t, a=p.b*v/p.t;
    cap+=s*exp(-s)*(em1(e)+exp(e)*sinhcm1(a));
    float u=2.*s/(cEdge*lambda*lambda);
    edge+=exp(-s)*em1((lambda-p.b)*u-.5*u*u);
 }
 float logTail=lp1(tail/normal);
 float logCap=lp1(64.*cap/float(cfg.panels));
 float margin=lp1(p.share*tail/normal+(1.-p.share)*64.*edge/float(cfg.panels))-logCap;
 out[idx]=float4(margin,logTail,logCap,lambda*margin-off*p.share*L/2.);
}
