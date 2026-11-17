/* Exact, gauge-reduced test for signed 1_S on S=supp(g) in F_2^7.
   The first eight affine-independent support signs are fixed +1.  Walsh at
   zero leaves only m=8,16,24 negative signs among the other 24 positions.
   Compile: clang -O3 -std=c11 signed_plateaued_enum.c -o signed_plateaued_enum */
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#define N 128
#define M 32
#define FIXED 8
#define FREE 24

static int parity7(uint8_t x) { return __builtin_parity((unsigned)x); }
static int gfun(uint8_t x) {
  int b0=(x>>0)&1,b1=(x>>1)&1,b2=(x>>2)&1,b3=(x>>3)&1;
  int b4=(x>>4)&1,b5=(x>>5)&1,b6=(x>>6)&1;
  return (b0*b1*b3) ^ (b0*b4*b5) ^ (b1*b2*b6);
}
static int rank7(const uint8_t *v, int n) {
  uint8_t piv[7]={0}; int r=0;
  for (int j=0;j<n;j++) { uint8_t x=v[j];
    for (int b=6;b>=0;b--) if ((x>>b)&1) {
      if (piv[b]) x^=piv[b]; else { piv[b]=x; r++; break; }
    }
  }
  return r;
}
static int acceptable(uint32_t negative, const uint32_t chi[N]) {
  for (int z=0;z<N;z++) {
    int w=M-2*__builtin_popcount(negative^chi[z]);
    if (w!=0 && w!=16 && w!=-16) return 0;
  }
  return 1;
}
static uint64_t choose(unsigned n,unsigned k) {
  uint64_t r=1; if(k>n-k)k=n-k;
  for(unsigned i=1;i<=k;i++) r=r*(n-k+i)/i;
  return r;
}
int main(void) {
  uint8_t support[M]; int ns=0;
  for(int x=0;x<N;x++) if(gfun((uint8_t)x)) support[ns++]=(uint8_t)x;
  if(ns!=M) { fprintf(stderr,"support size=%d, expected 32\n",ns); return 3; }
  int take[M]={0}; uint8_t basis[FIXED],d[7]; basis[0]=support[0]; take[0]=1; int nb=1,nd=0;
  for(int i=1;i<M && nb<FIXED;i++) {
    uint8_t tmp[7]; for(int j=0;j<nd;j++)tmp[j]=d[j]; tmp[nd]=support[i]^basis[0];
    if(rank7(tmp,nd+1)>nd) { d[nd]=tmp[nd]; nd++; basis[nb]=support[i]; nb++; take[i]=1; }
  }
  if(nb!=FIXED || rank7(d,7)!=7) { fprintf(stderr,"affine span rank failure\n"); return 4; }
  uint8_t p[M]; int at=0;
  for(int i=0;i<FIXED;i++)p[at++]=basis[i];
  for(int i=0;i<M;i++)if(!take[i])p[at++]=support[i];
  uint32_t chi[N]={0};
  for(int z=0;z<N;z++)for(int i=0;i<M;i++)if(parity7((uint8_t)z&p[i]))chi[z]|=UINT32_C(1)<<i;
  uint64_t checked=0, anticipated=2*choose(FREE,8)+1, solution_count=0; uint32_t witness=0;
  clock_t begin=clock();
  for(unsigned weight_index=0;weight_index<3;weight_index++) {
    unsigned m=(weight_index==0?8:(weight_index==1?16:24));
    uint32_t free_mask = (m==24)?((UINT32_C(1)<<FREE)-1):((UINT32_C(1)<<m)-1);
    for(;;) {
      uint32_t negative=free_mask<<FIXED; checked++;
      if(acceptable(negative,chi)) { if(solution_count==0) witness=negative; solution_count++; }
      if(m==24) break;
      uint32_t c=free_mask & -free_mask;
      uint32_t next=free_mask + c;
      if(next >= (UINT32_C(1)<<FREE)) break;
      free_mask=next | (((next^free_mask)>>2)/c);
    }
  }
  double sec=(double)(clock()-begin)/(double)CLOCKS_PER_SEC;
  printf("{\"schema\":\"a28-r032-signed-plateaued-enum-v2\",\"support_weight\":%d,\"affine_span_rank\":7,\"gauge_fixed_points\":[",ns);
  for(int i=0;i<FIXED;i++)printf("%s%u",i?",":"",(unsigned)basis[i]);
  printf("],\"states_target\":%llu,\"states_checked\":%llu,\"solution_count\":%llu,\"found\":%s,\"negative_sign_mask\":%u,\"seconds\":%.9f,\"support\":[",
    (unsigned long long)anticipated,(unsigned long long)checked,(unsigned long long)solution_count,solution_count?"true":"false",witness,sec);
  for(int i=0;i<M;i++)printf("%s%u",i?",":"",(unsigned)p[i]);
  printf("]}\n");
  return 0;
}
