#ifdef __METAL_VERSION__
#include <metal_stdlib>
using namespace metal;
typedef uint Word;
#define T thread
static Word bits(Word x){return popcount(x);}
static Word low_index(Word x){return ctz(x);}
#else
#include <stdint.h>
typedef uint32_t Word;
#define T
static Word bits(Word x){return __builtin_popcount(x);}
static Word low_index(Word x){return __builtin_ctz(x);}
#endif

struct Instance {
    Word a,b,family;
    Word table[16];
};
static Word draw(T Word *x){
    Word w=*x;
    w^=w<<13;w^=w>>17;w^=w<<5;
    *x=w;
    return w;
}
static Word maximum(Word a,Word b){return a>b?a:b;}
static Word minimum(Word a,Word b){return a<b?a:b;}
static Word sign_bit(Word x){return bits(x)&1u;}
static struct Instance make_instance(Word id){
    struct Instance m;
    m.family=id&3u;
    m.a=3+((id>>2)&1u);
    m.b=3+((id>>3)&1u);
    if(id==0){m.a=4;m.b=4;}
    if(id==1){m.a=3;m.b=3;}
    Word A=1u<<m.a,B=1u<<m.b;
    Word rng=id^0xa2806701u;
    if(!rng)rng=1;
    for(Word y=0;y<A;y++){
        Word type=m.family==3 ? draw(&rng)%3 : m.family;
        Word shift=draw(&rng)&(B-1),u=1,v=2,phase=draw(&rng)&1u;
        if(type==1){
            do{u=draw(&rng)&(B-1);}while(!u);
            do{v=draw(&rng)&(B-1);}while(!v||v==u);
        }
        if(id==0){shift=y;phase=0;type=0;}
        if(id==1){
            Word i=y%3;
            shift=0;phase=0;type=1;
            u=i==2?2:1;
            v=i==0?2:4;
        }
        Word table=0;
        if(type==2)table=draw(&rng)&((1u<<B)-1);
        else for(Word x=0;x<B;x++){
            Word z=sign_bit(shift&x)^phase;
            if(type==1)z^=sign_bit(u&x)&sign_bit(v&x);
            table|=z<<x;
        }
        m.table[y]=table;
    }
    return m;
}

static Word rank(Word points){
    Word span=1,r=0;
    while(points){
        Word v=low_index(points);
        points&=points-1;
        if((span>>v)&1u)continue;
        Word old=span;
        while(old){
            span|=1u<<(low_index(old)^v);
            old&=old-1;
        }
        r++;
    }
    return r;
}
static Word clique(Word vertices,T const Word *adj){
    Word stack[16]={vertices},depth=0,best=0;
    for(;;){
        Word active=stack[depth];
        if(depth+bits(active)<=best){
            if(depth==0)return best;
            depth--;
            continue;
        }
        Word node=low_index(active),rest=active&(active-1);
        stack[depth]=rest;
        depth++;
        best=maximum(best,depth);
        stack[depth]=rest&adj[node];
    }
}
static Word relate(Word labels,T Word *adj,Word check_all_pairs){
    Word count[16]={0},bad=0;
    for(Word p=labels;p;p&=p-1)for(Word q=p&(p-1);q;q&=q-1)count[low_index(p)^low_index(q)]++;
    for(Word p=labels;p;p&=p-1)for(Word q=p&(p-1);q;q&=q-1){
        Word x=low_index(p),y=low_index(q);
        if(count[x^y]<=1)bad|=check_all_pairs;
        else{
            adj[x]|=1u<<y;
            adj[y]|=1u<<x;
        }
    }
    return bad;
}

enum { RANK=0,TAIL=1,BASES=2,MAX_OCC_GAP=3,MAX_BASIS_GAP=4,ACT_CLIQUE=5,PROJ_CLIQUE=6,
       MIN_OCC=7,MIN_ACT=8,MIN_PROJ=9,FAMILY=10,SIG=11,A=12,B=13,
       GRAPH_HIST=14,BASIS_HIST=29,OCC_HIST=33,VALID=37,WIDTH=38 };
static void inspect(struct Instance m,T const Word *slices,T const Word *bands,T Word *rows){
    Word projected=0,check_union=0,actual[16]={0},global[16]={0},bad=0,q=0;
    Word ea=1u<<m.a,eb=1u<<m.b;
    for(Word y=0;y<ea;y++){
        check_union|=slices[y];
        bad|=relate(slices[y],actual,1);
    }
    for(Word beta=0;beta<eb;beta++)if(bands[beta]){
        projected|=1u<<beta;
        if(beta)q+=bands[beta];
    }
    bad|=projected!=check_union;
    relate(projected|1u,global,0);
    Word nonzero=projected&~1u,n=bits(nonzero),r=rank(nonzero);
    Word wa=clique(nonzero,actual),wp=clique(nonzero,global);
    if(wa>wp)return;
    for(Word i=0;i<eb;i++)bad|=(actual[i]&~global[i])!=0;
    rows[RANK]=r;
    rows[TAIL]=q;
    rows[ACT_CLIQUE]=wa;
    rows[PROJ_CLIQUE]=wp;
    rows[GRAPH_HIST+wp-wa]=1;
    rows[MIN_OCC]=rows[MIN_ACT]=rows[MIN_PROJ]=r;
    rows[FAMILY]=m.family;
    rows[A]=m.a;
    rows[B]=m.b;
    for(Word j=0;j<ea;j++)rows[SIG]=rows[SIG]*16777619u^(slices[j]^m.table[j]);
    for(Word j=0;j<eb;j++)rows[SIG]=rows[SIG]*16777619u^(bands[j]^(actual[j]<<8)^(global[j]<<16));
    if(r==0){
        rows[BASES]=rows[BASIS_HIST]=rows[OCC_HIST]=1;
        rows[VALID]=(bad==0);
        return;
    }
    Word labels[15],j=0;
    for(Word t=nonzero;t;t&=t-1)labels[j++]=low_index(t);
    for(Word combo=(1u<<r)-1;combo<(1u<<n);){
        Word mask=0,budget=0;
        for(Word t=combo;t;t&=t-1){
            Word beta=labels[low_index(t)];
            mask|=1u<<beta;
            budget+=bands[beta];
        }
        if(rank(mask)==r){
            Word occupancy=0;
            for(Word y=0;y<ea&&occupancy<r;y++)occupancy=maximum(occupancy,bits(mask&slices[y]));
            Word a=clique(mask,actual),p=clique(mask,global);
            if(occupancy>a || a>p || r*r>q*occupancy || r*r>budget*occupancy)return;
            rows[BASES]++;
            rows[MIN_OCC]=minimum(rows[MIN_OCC],occupancy);
            rows[MIN_ACT]=minimum(rows[MIN_ACT],a);
            rows[MIN_PROJ]=minimum(rows[MIN_PROJ],p);
            rows[MAX_OCC_GAP]=maximum(rows[MAX_OCC_GAP],a-occupancy);
            rows[MAX_BASIS_GAP]=maximum(rows[MAX_BASIS_GAP],p-a);
            rows[BASIS_HIST+p-a]++;
            rows[OCC_HIST+a-occupancy]++;
        }
        Word ripple=combo+(combo&(-combo));
        combo=ripple|((combo^ripple)>>(low_index(combo)+2));
    }
    rows[VALID]=(bad==0 && rows[BASES]>0);
}
