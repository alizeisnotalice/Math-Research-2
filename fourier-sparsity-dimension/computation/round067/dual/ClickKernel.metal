struct Configuration { uint first; uint length; };
kernel void audit(device uint *out [[buffer(0)]],constant Configuration& c [[buffer(1)]],
                  uint id [[thread_position_in_grid]]) {
    if(id>=c.length)return;
    Instance m=make_instance(c.first+id);
    uint ea=1u<<m.a,eb=1u<<m.b;
    uint active[16]={0},bands[16]={0},row[WIDTH]={0};
    int value[256];
    for(uint j=0;j<ea;j++){
        uint offset=j*eb;
        for(uint x=0;x<eb;x++)value[offset+x]=1-2*int((m.table[j]>>x)&1);
        for(uint d=1;d<eb;d*=2)for(uint x=0;x<eb;x++)if(!(x&d)){
            int a=value[offset+x],b=value[offset+x+d];
            value[offset+x]=a+b;
            value[offset+x+d]=a-b;
        }
        for(uint x=0;x<eb;x++)if(value[offset+x])active[j]|=1u<<x;
    }
    for(uint d=1;d<ea;d*=2)for(uint j=0;j<ea;j++)if(!(j&d))for(uint x=0;x<eb;x++){
        uint i=j*eb+x,k=(j+d)*eb+x;
        int a=value[i],b=value[k];
        value[i]=a+b;
        value[k]=a-b;
    }
    for(uint j=0;j<ea;j++)for(uint x=0;x<eb;x++)bands[x]+=(value[j*eb+x]!=0);
    inspect(m,active,bands,row);
    for(uint j=0;j<WIDTH;j++)out[id*WIDTH+j]=row[j];
}
