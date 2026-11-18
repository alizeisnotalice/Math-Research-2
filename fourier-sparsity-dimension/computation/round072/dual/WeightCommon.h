#ifdef __METAL_VERSION__
#include <metal_stdlib>
using namespace metal;
typedef uint Word;
typedef ulong Wide;
#define T thread
static Word bits(Word x){return popcount(x);}
static Word low_index(Word x){return ctz(x);}
#else
#include <stdint.h>
typedef uint32_t Word;
typedef uint64_t Wide;
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

static Word increment_rank(Word v,T Word *basis,Word width){
    for(Word j=width;j--;){
        if(!((v>>j)&1u))continue;
        if(!basis[j]){basis[j]=v;return 1;}
        v^=basis[j];
    }
    return 0;
}

enum { DIMENSION=0, L=1, M=2, S=3, PACKING_L=4, HOLES=5,
       INTRINSIC=6, VARYING=9, INVERSE_ACTIVE=12, PACK_IMPROVES=15,
       SIG=19, ADDR=20, DATA=21, FAMILY=22, VALID=23, WIDTH=24 };
enum { DEN=720720, RECIP_DEN=840 };

static void inspect(Word id,struct Instance m,T const Word *activity,T const int *full,T Word *row){
    Word A=1u<<m.a,B=1u<<m.b,d=0,at[16]={0},packing[16]={0},size[16]={0},
         p[8]={0},bad=0,hash=0;
    Word occupancy[4][16]={{0}},selected[4][16]={{0}};
    Wide factors[4]={0},scales[4]={1,DEN,(Wide)DEN*RECIP_DEN,DEN};
    for(Word t=0;t<A;t++)for(Word b=0;b<B;b++){
        if((activity[t]>>b)&1u)at[b]|=1u<<t;
        if(full[t*B+b]){
            size[b]|=1u<<t;
            d+=increment_rank(t*B+b,p,m.a+m.b);
        }
        hash=hash*16777619u^(Word)full[t*B+b];
    }
    row[DIMENSION]=d;
    row[ADDR]=m.a;
    row[DATA]=m.b;
    row[FAMILY]=m.family;
    Word rng=id^0xa2807201u;
    if(!rng)rng=1;
    for(Word b=0;b<B;b++){
        if((size[b]!=0)!=(at[b]!=0))return;
        if(!size[b])continue;
        Word star=low_index(size[b]),direction[4]={0},dim=0;
        for(Word s=size[b];s;s&=s-1)dim+=increment_rank(low_index(s)^star,direction,m.a);
        Word M_beta=1u<<dim,n=bits(size[b]),activity_n=bits(at[b]);
        row[HOLES]|=n<M_beta;
        row[M]+=M_beta;
        row[S]+=n;
        Word phase=low_index(at[b]);
        for(Word t=0;t<A;t++){
            Word off=t^phase,include=1;
            for(Word i=0;i<m.a;i++)include&=!sign_bit(off&direction[i]);
            if(include)packing[b]|=1u<<t;
        }
        bad|=bits(packing[b])*M_beta!=A || (packing[b]&~at[b])!=0;
        Word stretch=1+draw(&rng)%8;
        Word intr=M_beta*(DEN/n);
        Word w[4]={1,intr,intr*stretch,DEN/activity_n};
        factors[0]+=M_beta;
        factors[1]+=n;
        factors[2]+=n*(RECIP_DEN/stretch);
        factors[3]+=M_beta*activity_n;
        for(Word t=0;t<A;t++){
            for(Word k=0;k<4;k++){
                if((at[b]>>t)&1u)occupancy[k][t]+=w[k];
                if((packing[b]>>t)&1u)selected[k][t]+=w[k];
            }
        }
        hash=(hash*16777619u)^size[b];
        hash=(hash*16777619u)^packing[b];
        hash=(hash*16777619u)^stretch;
    }
    Wide base=0;
    for(Word k=0;k<4;k++){
        Word active=0,cosets=0;
        for(Word t=0;t<A;t++){
            active=maximum(active,occupancy[k][t]);
            cosets=maximum(cosets,selected[k][t]);
        }
        if(k==0){row[L]=active;row[PACKING_L]=cosets;base=active*factors[0];}
        Wide act_product=(Wide)active*factors[k],pack_product=(Wide)cosets*factors[k];
        bad|=pack_product>act_product || ((Wide)d*d)*scales[k]>7*pack_product;
        if(k){
            Wide scaled_base=base*scales[k];
            Word cmp=act_product<scaled_base?0:(act_product==scaled_base?1:2);
            row[INTRINSIC+3*(k-1)+cmp]=1;
        }
        row[PACK_IMPROVES+k]=(pack_product<act_product);
        hash=(hash*16777619u)^active;
        hash=(hash*16777619u)^cosets;
        hash=(hash*16777619u)^(Word)factors[k];
    }
    row[SIG]=hash;
    row[VALID]=(bad==0);
}
