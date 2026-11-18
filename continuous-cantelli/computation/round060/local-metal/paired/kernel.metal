struct Result { float exact, leading, A, B; };
kernel void matched(constant float2* nodes [[buffer(0)]],
                    device Result* output [[buffer(1)]],
                    constant uint& order [[buffer(2)]],
                    device const Point* input [[buffer(3)]],
                    uint j [[thread_position_in_grid]]) {
    const Point p=input[j];
    const float3 ab=integrate(p,nodes,order,false);
    const float paired=integrate(p,nodes,order,true).z;
    const float a=p.k*p.h/p.r;
    const float target=p.m+p.eps;
    const float scale=2.0f*sqrt(target*0.1591549430918953f)*exp(-0.5f*a*a*(1.0f-1.0f/target));
    output[j]=Result{scale*paired,scale*(ab.x-p.k*ab.y),ab.x,ab.y};
}
