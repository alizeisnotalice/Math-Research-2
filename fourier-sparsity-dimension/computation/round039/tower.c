// Enumerate all degree-2 and degree-3 monomials, followed by all affine terms.
// The lifted five-dimensional slice seed uses variables 0..3 and pair 4,5.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
typedef uint64_t Mask;
static Mask monomial(int variables) {
    Mask word=0;
    for(int x=0;x<64;x++) if((x&variables)==variables) word|=UINT64_C(1)<<x;
    return word;
}
static int bent(Mask word,const Mask chars[64]){
    for(int i=0;i<64;i++){
        int weight=__builtin_popcountll(word^chars[i]);
        if(weight!=28 && weight!=36)return 0;
    }return 1;
}
int main(int argc,char**argv){
    if(argc<4 || argc>8)return 2;
    int d=5+argc-3; int size=1<<d;
    unsigned start=atoi(argv[1]),step=atoi(argv[2]);
    Mask chars[64]={0},terms[35],phase[1024],gen[10];
    for(int i=0;i<64;i++)for(int x=0;x<64;x++)
        if(__builtin_parity(i&x))chars[i]|=UINT64_C(1)<<x;
    int n=0;
    for(int i=0;i<64;i++)if(__builtin_popcount(i)==3)terms[n++]=monomial(i);
    for(int i=0;n<35;i++)if(__builtin_popcount(i)==2)terms[n++]=monomial(i);
    Mask row=monomial(8),slice=~monomial(32);
    gen[0]=~row&slice;gen[1]=gen[0]&monomial(1);gen[2]=gen[0]&monomial(2);
    gen[3]=row&slice&monomial(1);gen[4]=row&slice&monomial(2);
    for(int i=3;i<argc;i++)gen[5+i-3]=strtoull(argv[i],NULL,16);
    phase[0]=monomial(3)^monomial(12)^monomial(48);
    for(int i=1;i<size;i++)phase[i]=phase[i&(i-1)]^gen[__builtin_ctz(i)];
    Mask support=0;for(int i=0;i<d;i++)support|=gen[i];
    unsigned long long kept=0, independent=0;int best=65;Mask witness=0;
    for(unsigned prefix=start;prefix<(1u<<18);prefix+=step){
        Mask w=0;for(int j=0;j<18;j++)if(prefix>>j&1)w^=terms[17+j];
        for(unsigned j=0;j<(1u<<17);j++){
            if(j)w^=terms[__builtin_ctz(j)];
            int i=0;while(i<size && bent(w^phase[i],chars))i++;
            if(i<size)continue;
            for(int linear=0;linear<64;linear++)for(int c=0;c<2;c++){
                Mask candidate=w^chars[linear]^(c?~UINT64_C(0):0);
                kept++;
                for(i=0;i<size;i++)if(candidate==(phase[i]^phase[0]))break;
                if(i==size){
                    independent++;
                    int size=__builtin_popcountll(candidate|support);
                    if(size<best){best=size;witness=candidate;}
                }
            }
        }
    }
    printf("%u %u %llu %llu %d %016llx\n",start,step,kept,independent,best,(unsigned long long)witness);
}
