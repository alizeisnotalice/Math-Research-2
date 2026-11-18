#ifdef __METAL_VERSION__
#include <metal_stdlib>
using namespace metal;
typedef uint Word;
#define T thread
#define G threadgroup
static Word bits(Word x){return popcount(x);}
static Word low(Word x){return ctz(x);}
#else
#include <stdint.h>
typedef uint32_t Word;
#define T
#define G
static Word bits(Word x){return __builtin_popcount(x);}
static Word low(Word x){return __builtin_ctz(x);}
#endif

struct Model{ Word v,y,family,labels[5],label_mask; };
struct Local{ int coeff[5]; Word zero_mode,phase; };

static Word draw(T Word *x){
    Word w=*x;w^=w<<13;w^=w>>17;w^=w<<5;*x=w;return w;
}
static Word sign_parity(Word x){return bits(x)&1u;}
static Word add_rank(Word v,T Word *p,Word width){
    for(Word j=width;j--;){
        if(!((v>>j)&1u))continue;
        if(!p[j]){p[j]=v;return 1;}
        v^=p[j];
    }
    return 0;
}
static struct Model make_model(Word id){
    struct Model m;
    m.v=2+(id>>2)%3;
    m.y=1+((id>>4)&1u);
    m.family=id&3u;
    if(!id){m.v=4;m.y=2;}
    Word s=id^0xa2806901u;
    if(!s)s=1;
    m.labels[0]=draw(&s)&15u;
    m.label_mask=1u<<m.labels[0];
    Word p[4]={0};
    for(Word j=1;j<5;j++){
        Word d;
        do{d=draw(&s)&15u;}while(!add_rank(d,p,4));
        m.labels[j]=m.labels[0]^d;
        m.label_mask|=1u<<m.labels[j];
    }
    return m;
}
static struct Local local_model(Word id,struct Model m,Word y,Word v){
    struct Local f={{0},0,0};
    Word s=id^0x691b156du^((y+1)*0x9e3779b9u)^((v+1)*0x846ca68bu);
    if(!s)s=1;
    Word omitted=draw(&s)%5,amount=0;
    if(m.family==0){
        for(Word i=0;i<5;i++)if(i!=omitted)f.coeff[i]=1;
    }else{
        amount=(m.family==2)?2:2*(draw(&s)%3);
        for(Word k=0;k<amount;k++)f.coeff[draw(&s)%5]++;
    }
    Word phases=draw(&s);
    for(Word i=0;i<5;i++)if((phases>>i)&1u)f.coeff[i]=-f.coeff[i];
    f.zero_mode=(m.family!=0 && m.family!=1 && ((phases>>7)&1u));
    f.phase=(phases>>8)&1u;
    return f;
}
static int control_mean(struct Model m,struct Local f,Word u){
    int g=0;
    for(Word i=0;i<5;i++)g+=sign_parity(m.labels[i]&u)?-f.coeff[i]:f.coeff[i];
    return g;
}
static Word negative_mask(int g,struct Local f){
    if(g==4)return 0;
    if(g==-4)return 0xffffffffu;
    if(g==2)return 0x88888888u;
    if(g==-2)return 0x77777777u;
    Word value=f.zero_mode?0xffff0000u:0x0fffF000u;
    return f.phase?~value:value;
}

static void transform(T int *a,Word len){
    for(Word stride=1;stride<len;stride*=2)
        for(Word off=0;off<len;off+=2*stride)
            for(Word i=0;i<stride;i++){
                int l=a[off+i],r=a[off+i+stride];
                a[off+i]=l+r;a[off+i+stride]=l-r;
            }
}

enum { BATCH=0,VDIM=1,YDIM=2,COMMON=3,AFFINE_RANK=4,EXCESS=5,FIBER_MAX=6,SHARP=7,ACTUAL_SHARP=8,
       THIST=9,ENV_MAX=14,ACTUAL_MAX=15,ZERO=16,HALF=17,ONE=18,FINGERPRINT=19,FAMILY=20,VALID=21,WIDTH=22 };

static Word private_cost(int measured,Word support,int prescribed,T Word *bad){
    Word count=bits(support),mag=measured<0?-measured:measured;
    Word haszero=support&1u;
    Word nonzero=support&~1u,p[5]={0},r=0;
    for(Word set=nonzero;set;set&=set-1)r+=add_rank(low(set),p,5);
    *bad|=(mag!=0 && mag!=2 && mag!=4) || measured!=prescribed;
    *bad|=((mag==4) && (support!=1u)) ||
          ((mag==2) && (count!=4||!haszero||r!=2)) ||
          ((mag==0) && (haszero||!((count==1&&r==1)||(count==4&&r==3))));
    Word n=bits(nonzero);
    if(!n)return 0;
    *bad|=(12*r*r)%n!=0;
    return 12*r*r/n;
}
static Word validate_bag(Word mask){
    Word pairs[32]={0};
    for(Word p=mask;p;p&=p-1)for(Word q=p&(p-1);q;q&=q-1)pairs[low(p)^low(q)]++;
    for(Word p=mask;p;p&=p-1)for(Word q=p&(p-1);q;q&=q-1)for(Word r=q&(q-1);r;r&=r-1){
        Word a=low(p),b=low(q),c=low(r),d=a^b^c;
        if(d>c && ((mask>>d)&1u) && pairs[a^b]!=2 && pairs[a^c]!=2 && pairs[a^d]!=2)return 1;
    }
    return 0;
}

static Word inspect_address(struct Model m,G const int *means,G const Word *types,G const Word *support,
                            T Word *core,T Word *row){
    Word V=1u<<m.v,bad=0,bag=0;
    int urows[256];
    for(Word v=0;v<V;v++){
        Word env=0,actual=0,t=0,class_count[3]={0};
        T int *coeff=urows+(v<<4);
        for(Word u=0;u<16;u++){
            Word point=(v<<4)|u;
            int g=means[point];
            Word mag=g<0?-g:g;
            if(mag!=0&&mag!=2&&mag!=4)return 1;
            Word upper=mag==4?0:(mag==2?16:27);
            coeff[u]=g;
            class_count[mag/2]++;
            env+=upper;actual+=types[point];
            bad|=types[point]>upper;
            bag|=support[point];
            row[FINGERPRINT]=row[FINGERPRINT]*16777619u^((Word)g^support[point]);
        }
        transform(coeff,16);
        for(Word beta=0;beta<16;beta++)if(coeff[beta]){
            t++;
            bad|=!((m.label_mask>>beta)&1u) || coeff[beta]%16!=0;
        }
        if(t>4)return 1;
        Word e=192*t+env,a=192*t+actual;
        bad|=(e>1058||a>e);
        if(t==4)bad|=(class_count[0]!=6||class_count[1]!=8||class_count[2]!=2||e!=1058);
        row[BATCH]++;
        row[THIST+t]++;
        row[SHARP]+=(e==1058);
        row[ACTUAL_SHARP]+=(a==1058);
        if(e>row[ENV_MAX])row[ENV_MAX]=e;
        if(a>row[ACTUAL_MAX])row[ACTUAL_MAX]=a;
        row[ZERO]+=class_count[0];
        row[HALF]+=class_count[1];
        row[ONE]+=class_count[2];
    }
    bad|=validate_bag(bag|1u);
    for(Word b=0;b<16;b++){
        int line[16];
        for(Word v=0;v<V;v++)line[v]=urows[(v<<4)|b];
        transform(line,V);
        for(Word v=0;v<V;v++)if(line[v])core[b]|=1u<<v;
    }
    return bad;
}
static Word finish(struct Model m,T const Word *core,T Word *row){
    Word bad=0,anchor=0xffffu,p[8]={0},r=0,c=0;
    for(Word beta=0;beta<16;beta++){
        bad|=core[beta] && !((m.label_mask>>beta)&1u);
        Word size=bits(core[beta]);
        if(size>row[FIBER_MAX])row[FIBER_MAX]=size;
        c+=size;
        for(Word v=core[beta];v;v&=v-1){
            Word frequency=(low(v)<<4)|beta;
            if(anchor==0xffffu)anchor=frequency;
            r+=add_rank(frequency^anchor,p,m.v+4);
        }
        row[FINGERPRINT]=row[FINGERPRINT]*16777619u^core[beta];
    }
    row[VDIM]=m.v;row[YDIM]=m.y;row[FAMILY]=m.family;
    row[COMMON]=c;row[AFFINE_RANK]=r;
    row[EXCESS]=c ? c-1-r:0;
    return bad;
}
