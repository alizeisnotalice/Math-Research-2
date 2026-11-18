#ifdef __METAL_VERSION__
#include <metal_stdlib>
using namespace metal;
typedef uint Word;
typedef ulong Wide;
#define T thread
static Word parity(Word x){return popcount(x)&1u;}
#else
#include <stdint.h>
typedef uint32_t Word;
typedef uint64_t Wide;
#define T
static Word parity(Word x){return __builtin_parity(x);}
#endif

struct Model { Word k,n,family; };
struct Leaf { Word n,forms[4],shift,phase,table; };
static Word draw(T Word *x){
    Word w=*x;w^=w<<13;w^=w>>17;w^=w<<5;*x=w;return w;
}
static Word add_rank(Word v,T Word *p,Word size){
    for(Word j=size;j--;){
        if(!((v>>j)&1u))continue;
        if(!p[j]){p[j]=v;return 1;}
        v^=p[j];
    }
    return 0;
}
static Word remainder(Word v,T const Word *p,Word size){
    for(Word j=size;j--;)if((v>>j)&1u)v^=p[j];
    return v;
}
static struct Model make_model(Word id){
    struct Model m;
    m.k=1+((id>>2)&1u);
    m.n=5+((id>>3)&1u);
    m.family=id&3u;
    if(id<6){m.n=5;m.k=id==0?1:2;}
    return m;
}
static struct Leaf leaf(Word id,struct Model m,Word address){
    struct Leaf f={0,{0},0,0,0};
    if(id<6){
        Word chars[6][4]={{0,0,0,0},{4,4,1,5},{1,2,3,4},{0,1,2,4},{7,7,2,2},{1,1,1,2}};
        f.shift=chars[id][address];
        f.phase=(id<4)? (address&1u):0;
        if(id==0){
            f.n=3;f.table=1;
            for(Word j=0;j<3;j++)f.forms[j]=1u<<(j+2*address);
        }
        if(id==3&&address==0){f.n=2;f.table=8;f.forms[0]=1;f.forms[1]=2;}
        return f;
    }
    Word seed=id^0xa2807401u^((address+1)*0x846ca68bu);
    if(!seed)seed=1;
    Word type=m.family==3?draw(&seed)%4:m.family,bits=(1u<<m.n)-1,p[6]={0};
    f.shift=draw(&seed)&bits;
    f.phase=draw(&seed)&1u;
    if(type==0)return f;
    f.n=type==1?2:3+(draw(&seed)&1u);
    f.table=type==1?8:(type==2?1:(draw(&seed)&((1u<<(1u<<f.n))-1)));
    for(Word j=0;j<f.n;j++){
        Word d;
        do {d=draw(&seed)&bits;} while(!add_rank(d,p,m.n));
        f.forms[j]=d;
    }
    return f;
}
static Word evaluate(struct Leaf f,Word x){
    Word target=0;
    for(Word i=0;i<f.n;i++)target|=parity(x&f.forms[i])<<i;
    return ((f.table>>target)^f.phase^parity(f.shift&x))&1u;
}

static void transform(T int *row,Word count){
    for(Word stride=1;stride<count;stride*=2)
        for(Word k=0;k<count;k+=2*stride)
            for(Word j=0;j<stride;j++){
                int l=row[k+j],r=row[k+j+stride];
                row[k+j]=l+r;row[k+j+stride]=l-r;
            }
}

enum { K=0,DATA=1,PARENT_RANK=2,DATA_RANK=3,PARENT_CORE_RANK=4,CORE_RANK=5,
       PARENT_POINTS=6,TAIL_POINTS=7,BANDS=8,REPEATED=9,EXCLUSIVE=10,INNER_EXCLUSIVE=11,
       TAIL_BANDS=12,CHILD_RANK_SUM=13,TAIL_RANK=14,GLUE_GAP=15,EQUAL_BANDS=16,STRICT_BANDS=17,
       NO_TAIL=18,ADDRESS_DEFICIT=19,SIG=20,FAMILY=21,GLUE_STRICT=22,HIST=23,VALID=39,WIDTH=40 };

static void inspect(struct Model m,T const int *child,T const int *parent,T Word *row){
    Word B=1u<<m.k,Z=1u<<m.n;
    Word core[6]={0},all_beta[6]={0},core_rank=0,r_data=0,r_parent=0,r_core_parent=0,r_tail=0;
    Word p[8]={0},pc[8]={0},tail[6]={0},perchild[4][6]={{0}},glued[8]={0},r_glued=0;
    Word participants[64]={0},sum_s=0,s_tail=0;
    for(Word beta=0;beta<Z;beta++){
        for(Word b=0;b<B;b++)if(child[b*Z+beta])participants[beta]|=1u<<b;
        Word active=participants[beta];
        if(active){
            r_data+=add_rank(beta,all_beta,m.n);
            if(active&(active-1))core_rank+=add_rank(beta,core,m.n);
        }
    }
    row[K]=m.k;row[DATA]=m.n;row[FAMILY]=m.family;
    Word hash=0;
    for(Word beta=0;beta<Z;beta++){
        Word active=participants[beta],a=0,m_band=0,rep=B;
        for(Word b=0;b<B;b++)a+=(active>>b)&1u;
        Word reduced=remainder(beta,core,m.n),outside=reduced!=0;
        if(outside&&a>1)return;
        if(outside&&a==1){
            row[TAIL_BANDS]++;
            r_tail+=add_rank(reduced,tail,m.n);
            for(Word b=0;b<B;b++)if((active>>b)&1u)row[CHILD_RANK_SUM]+=add_rank(reduced,perchild[b],m.n);
        }
        for(Word alpha=0;alpha<B;alpha++){
            int expect=0;
            for(Word b=0;b<B;b++){
                int w=child[b*Z+beta];
                expect+=parity(alpha&b)?-w:w;
            }
            int w=parent[alpha*Z+beta];
            if(w!=expect)return;
            hash=hash*16777619u^(Word)w;
            if(w){
                m_band++;
                r_parent+=add_rank(alpha*Z+beta,p,m.k+m.n);
                if(!outside)r_core_parent+=add_rank(alpha*Z+beta,pc,m.k+m.n);
                if(rep==B){
                    rep=alpha;
                    r_glued+=add_rank(alpha*Z+beta,glued,m.k+m.n);
                }else r_glued+=add_rank((alpha^rep)*Z,glued,m.k+m.n);
            }
        }
        if(!a){if(m_band)return;continue;}
        sum_s+=m_band;
        if(m_band*a<B)return;
        row[HIST+4*(a-1)+m_band-1]++;
        row[EQUAL_BANDS]+=(m_band*a==B);
        row[STRICT_BANDS]+=(m_band*a>B);
        row[BANDS]++;
        if(a==1){
            if(m_band!=B)return;
            row[EXCLUSIVE]++;
            row[INNER_EXCLUSIVE]+=!outside;
            if(outside)s_tail+=m_band;
        }else row[REPEATED]++;
    }
    if(r_parent!=r_glued || core_rank+r_tail!=r_data || r_parent<r_data ||
       r_parent>r_data+m.k || r_core_parent<core_rank || r_core_parent>core_rank+m.k ||
       r_tail>row[CHILD_RANK_SUM] || s_tail!=B*row[TAIL_BANDS] ||
       row[TAIL_BANDS]+row[INNER_EXCLUSIVE]!=row[EXCLUSIVE])return;
    if(row[EXCLUSIVE] && r_parent!=r_data+m.k)return;
    if(row[INNER_EXCLUSIVE] && r_core_parent!=core_rank+m.k)return;
    if(s_tail){
        if(r_parent!=r_core_parent+r_tail+m.k-(r_core_parent-core_rank))return;
    }else if(r_parent!=r_core_parent)return;
    row[PARENT_RANK]=r_parent;
    row[DATA_RANK]=r_data;
    row[PARENT_CORE_RANK]=r_core_parent;
    row[CORE_RANK]=core_rank;
    row[PARENT_POINTS]=sum_s;
    row[TAIL_POINTS]=s_tail;
    row[TAIL_RANK]=r_tail;
    row[GLUE_GAP]=row[CHILD_RANK_SUM]-r_tail;
    row[GLUE_STRICT]=row[GLUE_GAP]>0;
    row[NO_TAIL]=!s_tail;
    row[ADDRESS_DEFICIT]=r_data+m.k-r_parent;
    row[SIG]=hash;
    row[VALID]=1;
}
