struct Launch { uint first; };
constant uint LANES=128;

static void transform(threadgroup int *data,uint len,uint lane) {
    for(uint gap=1;gap<len;gap<<=1){
        for(uint p=lane;p<len/2;p+=LANES){
            uint a=(p/gap)*2*gap+p%gap,b=a+gap;
            int u=data[a],v=data[b];
            data[a]=u+v;data[b]=u-v;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
}

kernel void restriction_audit(device uint *records [[buffer(0)]],
                              constant Launch& launch [[buffer(1)]],
                              uint gid [[threadgroup_position_in_grid]],
                              uint lane [[thread_index_in_threadgroup]]) {
    Model model=create_model(launch.first+gid);
    uint M=1u<<model.h,V=1u<<model.k,N=8*M,small=8*V,scale=M/V;
    threadgroup int full[512],actual[256];
    threadgroup uint failed[LANES], projected[64];
    for(uint x=lane;x<N;x+=LANES)full[x]=evaluate(x,model);
    for(uint x=lane;x<small;x+=LANES)actual[x]=evaluate(restriction_point(x,model),model);
    if(lane<M)projected[lane]=image(lane,model);
    threadgroup_barrier(mem_flags::mem_threadgroup);
    transform(full,N,lane);
    transform(actual,small,lane);
    uint invalid=0;
    for(uint j=lane;j<small;j+=LANES){
        uint beta=j/V;
        int folded=0;
        for(uint a=0;a<M;a++)if(projected[a]==j%V){
            int c=full[beta*M+a];
            folded+=parity(a&model.shift)?-c:c;
        }
        invalid|=(folded!=actual[j]*int(scale));
    }
    failed[lane]=invalid;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if(lane!=0)return;

    uint result[WIDTH]={0},old_rank[3]={0},new_rank[3]={0};
    result[H]=model.h;
    result[K]=model.k;
    result[HOLED]=model.holes;
    for(uint i=0;i<LANES;i++)invalid|=failed[i];
    uint step=1u<<(model.n-model.h+1+model.holes);
    for(uint beta=0;beta<8;beta++){
        uint start=M,old_count=0,new_count=0;
        for(uint input=0;input<M;input++){
            int value=full[beta*M+input];
            uint mag=uint(abs(value));
            uint odd=(mag/step)&1u;
            uint pattern=core_pattern_weight(input,model);
            uint expected_odd=beta==0 && (!model.holes || (pattern&1u));
            invalid|=(mag%step!=0 || odd!=expected_odd);
            result[ODD_CORE]+=odd;
            if(!beta){
                invalid|=((value!=0)!=(pattern!=2));
                result[CORE]+=(value!=0);
            }
            if(value){
                if(start==M)start=input;
                old_count++;
            }
        }
        for(uint a=0;a<V;a++){
            int c=actual[beta*V+a];
            result[SPECTRAL_HASH]=result[SPECTRAL_HASH]*16777619u^uint(c);
            new_count+=(c!=0);
        }
        result[ORIG]+=old_count;
        result[RESTRICTED]+=new_count;
        if(beta==0)continue;
        if(new_count)result[D_AFTER]+=add_rank(beta,new_rank,3);
        if(!old_count||!add_rank(beta,old_rank,3))continue;
        result[D_BEFORE]++;
        result[MISSING]+=(new_count==0);
        result[ELIGIBLE]+=old_count-1;
        uint pairs=0;
        for(uint a=start+1;a<M;a++)if(full[beta*M+a]&&projected[a]==projected[start])pairs++;
        result[PAIRS]+=pairs;
        result[COLLIDING]+=(pairs!=0);
        invalid|=(!new_count&&!pairs);
    }
    invalid|=result[CORE]!=(model.holes?5*M/8:M) || result[ODD_CORE]!=(model.holes?M/2:M) || result[D_AFTER]>result[D_BEFORE] ||
        result[D_BEFORE]-result[D_AFTER]>result[MISSING] ||
        result[MISSING]>result[COLLIDING] || result[COLLIDING]>result[PAIRS] ||
        result[PAIRS]>result[ELIGIBLE] || result[ELIGIBLE]>result[ORIG]-result[CORE]-result[D_BEFORE];
    result[VALID]=(invalid==0);
    for(uint j=0;j<WIDTH;j++)records[gid*WIDTH+j]=result[j];
}
