struct Configuration { uint first,length; };
kernel void audit(device uint*out [[buffer(0)]],constant Configuration&c [[buffer(1)]],
                  uint id [[thread_position_in_grid]]){
    if(id>=c.length)return;
    Instance m=make_instance(c.first+id);
    uint A=1u<<m.a,B=1u<<m.b,live[16]={0},row[WIDTH]={0};
    int full[256];
    for(uint t=0;t<A;t++){
        uint offset=t*B;
        for(uint z=0;z<B;z++)full[offset+z]=1-2*int((m.table[t]>>z)&1);
        for(uint d=1;d<B;d*=2)for(uint z=0;z<B;z++)if(!(z&d)){
            int a=full[offset+z],b=full[offset+z+d];
            full[offset+z]=a+b;full[offset+z+d]=a-b;
        }
        for(uint z=0;z<B;z++)if(full[offset+z])live[t]|=1u<<z;
    }
    for(uint d=1;d<A;d*=2)for(uint t=0;t<A;t++)if(!(t&d))for(uint b=0;b<B;b++){
        uint i=t*B+b,j=(t+d)*B+b;
        int a=full[i],c=full[j];full[i]=a+c;full[j]=a-c;
    }
    inspect(c.first+id,m,live,full,row);
    for(uint i=0;i<WIDTH;i++)out[id*WIDTH+i]=row[i];
}
