struct Configuration { uint first,length; };
kernel void audit(device uint*out [[buffer(0)]],constant Configuration &c [[buffer(1)]],
                  uint index [[thread_position_in_grid]]){
    if(index>=c.length)return;
    uint id=c.first+index;
    Model m=make_model(id);
    uint B=1u<<m.k,Z=1u<<m.n,row[WIDTH]={0};
    int local[256],combined[256];
    for(uint b=0;b<B;b++){
        Leaf f=leaf(id,m,b);
        uint offset=b*Z;
        int restriction[64];
        for(uint x=0;x<Z;x++)restriction[x]=combined[offset+x]=1-2*int(evaluate(f,x));
        transform(restriction,Z);
        for(uint x=0;x<Z;x++)local[offset+x]=restriction[x];
    }
    transform(combined,B*Z);
    inspect(m,local,combined,row);
    for(uint j=0;j<WIDTH;j++)out[index*WIDTH+j]=row[j];
}
