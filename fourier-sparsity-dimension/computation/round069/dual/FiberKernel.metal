struct Configuration { uint first; };
kernel void audit(device uint *out [[buffer(0)]],constant Configuration& conf [[buffer(1)]],
                  uint id [[threadgroup_position_in_grid]],uint lane [[thread_index_in_threadgroup]]) {
    Model m=make_model(id+conf.first);
    uint V=1u<<m.v,Y=1u<<m.y;
    threadgroup Local sites[16];
    threadgroup uint supports[256],costs[256],invalid[256];
    threadgroup int means[256];
    uint core[16]={0},rows[WIDTH]={0},bad=0;
    for(uint y=0;y<Y;y++){
        if(lane<V)sites[lane]=local_model(id+conf.first,m,y,lane);
        threadgroup_barrier(mem_flags::mem_threadgroup);
        for(uint site=lane;site<V*16;site+=128){
            Local f=sites[site>>4];
            int expected=control_mean(m,f,site&15u);
            uint mask=negative_mask(expected,f),support=0;
            int z[32];
            for(uint x=0;x<32;x++)z[x]=1-2*int((mask>>x)&1u);
            transform(z,32);
            for(uint b=0;b<32;b++)if(z[b])support|=1u<<b;
            uint error=(z[0]%8!=0);
            means[site]=z[0]/8;
            costs[site]=private_cost(means[site],support,expected,&error);
            supports[site]=support;
            invalid[site]=error;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
        if(lane==0){
            for(uint site=0;site<V*16;site++)bad|=invalid[site];
            bad|=inspect_address(m,means,costs,supports,core,rows);
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if(lane==0){
        bad|=finish(m,core,rows);
        rows[VALID]=(bad==0);
        for(uint j=0;j<WIDTH;j++)out[id*WIDTH+j]=rows[j];
    }
}
