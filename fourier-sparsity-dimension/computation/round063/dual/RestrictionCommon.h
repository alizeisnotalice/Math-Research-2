#ifdef __METAL_VERSION__
#include <metal_stdlib>
using namespace metal;
typedef uint Word;
#define THREAD thread
#else
#include <stdint.h>
typedef uint32_t Word;
#define THREAD
#endif

struct Model {
    Word h, k, n, shift, vectors[5], forms[3], table;
};
static Word parity(Word x) {
#ifdef __METAL_VERSION__
    return popcount(x)&1;
#else
    return __builtin_parity(x);
#endif
}
static Word random_word(THREAD Word *state) {
    Word w=*state;
    w^=w<<13; w^=w>>17; w^=w<<5;
    *state=w;
    return w;
}
static Word add_rank(Word v, THREAD Word *basis, Word dimension) {
    for(Word j=dimension; j--;) {
        if(!((v>>j)&1))continue;
        if(!basis[j]){basis[j]=v;return 1;}
        v^=basis[j];
    }
    return 0;
}
static struct Model create_model(Word id) {
    struct Model m;
    m.h=5+((id>>1)&1);
    m.k=1+(id>>2)%(m.h-1);
    m.n=m.h+3;
    Word state=id^0xa2806301u;
    if(!state)state=1;
    Word mask=(1u<<m.h)-1;
    Word shape=mask^1u;
    m.shift=random_word(&state)&mask;
    m.forms[0]=random_word(&state)&shape;
    m.forms[1]=(random_word(&state)&shape)|(1u<<m.h);
    m.forms[2]=(random_word(&state)&shape)|(2u<<m.h);
    m.table=random_word(&state)&255u;
    Word b[6]={0};
    for(Word j=0;j<m.k;j++){
        Word v;
        do {v=random_word(&state)&mask;}while(!add_rank(v,b,m.h));
        m.vectors[j]=v;
    }
    if(id<4){
        m.shift=id;
        m.k=2;
        m.vectors[0]=2;
        m.vectors[1]=4;
        m.forms[0]=2;
        m.forms[1]=4|(1u<<m.h);
        m.forms[2]=8|(2u<<m.h);
        m.table=0x97;
    }
    return m;
}
static Word image(Word a,struct Model m) {
    Word out=0;
    for(Word j=0;j<m.k;j++)out|=parity(a&m.vectors[j])<<j;
    return out;
}
static Word restriction_point(Word index,struct Model m) {
    Word x=m.shift|((index>>m.k)<<m.h);
    for(Word j=0;j<m.k;j++)if((index>>j)&1)x^=m.vectors[j];
    return x;
}
static int evaluate(Word x,struct Model m) {
    if(!(x&1))return (x&((1u<<m.h)-2u)) ? 1 : -1;
    Word a=parity(x&m.forms[0])|(parity(x&m.forms[1])<<1)|(parity(x&m.forms[2])<<2);
    return ((m.table>>a)^ (x>>(m.n-1)))&1u ? -1 : 1;
}

enum { WIDTH=13, 
       ORIG=0, RESTRICTED=1, D_BEFORE=2, D_AFTER=3, ELIGIBLE=4, PAIRS=5, COLLIDING=6,
       MISSING=7, CORE=8, H=9, K=10, SPECTRAL_HASH=11, VALID=12 };
