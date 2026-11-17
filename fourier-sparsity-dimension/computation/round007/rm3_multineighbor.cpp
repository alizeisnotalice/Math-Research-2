// Exact one-hyperplane-neighbour search around corrected-results row 1.
// Search scope: one fixed W=C+7 relative generators (relative generators 4,5,6,7 kept), modulo W inside RM(3,6).  This is not a global A28 search.
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "neighbor_input.hpp"

using Mask = std::uint64_t;
constexpr int N=6, P=64;
constexpr Mask FULL=~Mask(0);
constexpr int KEEP_FROM=4;

static int pc(Mask x) { return __builtin_popcountll(x); }
static Mask phase_q() { Mask q=0; for(int x=0;x<P;x++) {
  const int bit=(pc((x&7)&INPUT_PERM[x>>3]) ^ ((INPUT_H_BITS>>(x>>3))&1))&1;
  if(bit) q|=Mask(1)<<x;
} return q; }
struct Space {
  std::array<Mask,64> piv{};
  int rank=0;
  bool insert(Mask x) { for(int b=63;b>=0 && x;b--) if((x>>b)&1) {
    if(piv[b]) x^=piv[b]; else {piv[b]=x;rank++;return true;}
  } return false; }
  bool contains(Mask x) const { for(int b=63;b>=0 && x;b--) if((x>>b)&1) {if(!piv[b]) return false; x^=piv[b];} return true; }
};
static Mask coordinate_mask(int j) { Mask r=0; for(int x=0;x<P;x++) if((x>>j)&1) r|=Mask(1)<<x; return r; }
static Mask monomial_mask(int variables) { Mask r=0; for(int x=0;x<P;x++) if((x&variables)==variables)r|=Mask(1)<<x; return r; }
static bool bent(Mask f) {
  int a[P]; for(int x=0;x<P;x++) a[x]=((f>>x)&1)?-1:1;
  for(int step=1;step<P;step<<=1) for(int s=0;s<P;s+=2*step) for(int j=s;j<s+step;j++) {const int x=a[j],y=a[j+step];a[j]=x+y;a[j+step]=x-y;}
  for(int x:a) if(x!=8 && x!=-8) return false;
  return true;
}
static Mask phase_for_index(std::uint64_t index,const std::vector<Mask>& qbasis,Mask base) {
  const std::uint64_t g=index^(index>>1); Mask out=base;
  for(int j=0;j<(int)qbasis.size();j++) if((g>>j)&1) out^=qbasis[j];
  return out;
}
static bool all_w_slices_bent(Mask f,const std::array<Mask,4>& relw) {
  Mask cur=f, prev=0;
  for(unsigned z=0;z<16;z++) {const unsigned gray=z^(z>>1); if(z) cur^=relw[__builtin_ctz(gray^prev)]; if(!bent(cur))return false; prev=gray;}
  return true;
}
static int h_in_span(const Space& w,Mask f) { Space s=w;s.insert(f);int n=0;for(Mask h:INPUT_H)if(s.contains(h))n++;return n; }
int main(int argc,char**argv) {
  if(argc!=4){std::cerr<<"usage: rm3_neighbor_search START COUNT PREFIX\n";return 2;}
  const std::uint64_t start=std::stoull(argv[1]), count=std::stoull(argv[2]); const std::string prefix=argv[3];
  if(start>(1ULL<<31) || count>(1ULL<<31)-start){std::cerr<<"range error\n";return 8;}
  Space w; w.insert(FULL); for(int j=0;j<N;j++)w.insert(coordinate_mask(j));
  std::array<Mask,4> relw{}; int at=0; for(int j=0;j<8;j++)if(j>=KEEP_FROM){w.insert(INPUT_REL[j]);relw[at++]=INPUT_REL[j];}
  if(w.rank!=11 || h_in_span(w,0)==8){std::cerr<<"invalid W\n";return 3;}
  Space rm; std::vector<Mask> monomials; for(int m=0;m<(1<<N);m++)if(pc(m)<=3){Mask z=monomial_mask(m);rm.insert(z);monomials.push_back(z);} if(rm.rank!=42){std::cerr<<"RM rank\n";return 4;}
  for(int j=0;j<64;j++)if(w.piv[j]&&!rm.contains(w.piv[j])){std::cerr<<"W not in RM3\n";return 5;}
  Space extend=w; std::vector<Mask> qbasis; for(Mask z:monomials)if(extend.insert(z))qbasis.push_back(z);
  if(qbasis.size()!=31 || extend.rank!=42){std::cerr<<"quotient rank\n";return 6;}
  const Mask f0=phase_q(); if(!rm.contains(f0) || !all_w_slices_bent(f0,relw)){std::cerr<<"baseline W illegal\n";return 7;}
  const auto began=std::chrono::steady_clock::now(); std::ofstream found(prefix+".found.jsonl");
  std::uint64_t hits=0; int max_h=0;
  Mask f=phase_for_index(start,qbasis,f0);
  for(std::uint64_t i=start;i<start+count;i++) { if(all_w_slices_bent(f,relw)) {const int hn=h_in_span(w,f^f0);max_h=std::max(max_h,hn);hits++;found<<"{\"index\":"<<i<<",\"phase_hex\":\"0x"<<std::hex<<f<<"\",\"relative_phase_hex\":\"0x"<<(f^f0)<<std::dec<<"\",\"H_blocks_in_W_plus_relative\":"<<hn<<",\"all_H_blocks\":"<<(hn==8?"true":"false")<<"}\n";} if(i+1<start+count) f^=qbasis[__builtin_ctzll(i+1)]; }
  const double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
  std::ofstream out(prefix+".summary.json");
  out<<"{\"status\":\"passed\",\"scope\":\"fixed_W_keep_relative_generators_4_5_6_7_RM3_mod_W\",\"start\":"<<start<<",\"count\":"<<count<<",\"tested\":"<<count<<",\"hits\":"<<hits<<",\"max_H_blocks_in_hit\":"<<max_h<<",\"seconds\":"<<sec<<",\"W_rank\":11,\"RM3_rank\":42,\"quotient_rank\":31,\"kept_relative_indices\":[4,5,6,7],\"H_blocks_in_W\":"<<h_in_span(w,0)<<"}\n";
  std::cout<<"status=passed start="<<start<<" count="<<count<<" hits="<<hits<<" seconds="<<sec<<" H_in_W="<<h_in_span(w,0)<<"\n";
}
