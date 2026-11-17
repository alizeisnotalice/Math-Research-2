// Independent exact enumeration: reverse monomial complement and autocorrelation.
#include <array>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "neighbor_input.hpp"
using U=std::uint64_t;
struct Basis {
 std::array<U,64> p{};int dim=0;
 U reduce(U x)const{for(int j=63;j>=0;--j)if((x>>j)&1)x^=p[j];return x;}
 bool add(U x){x=reduce(x);if(!x)return false;int j=63-__builtin_clzll(x);p[j]=x;++dim;return true;}
};
U monomial(unsigned mask){U r=0;for(unsigned x=0;x<64;++x)if((x&mask)==mask)r|=U(1)<<x;return r;}
std::array<U,6> low{};
U shift(U f,unsigned a){for(int j=0;j<6;++j)if((a>>j)&1){unsigned s=1U<<j;f=((f&low[j])<<s)|((f&~low[j])>>s);}return f;}
bool bent_derivative(U f){for(unsigned a=1;a<64;++a)if(__builtin_popcountll(f^shift(f,a))!=32)return false;return true;}
bool bent_walsh(U f){int a[64];for(int i=0;i<64;++i)a[i]=((f>>i)&1)?-1:1;for(int s=1;s<64;s*=2)for(int b=0;b<64;b+=s*2)for(int j=b;j<b+s;++j){int x=a[j],y=a[j+s];a[j]=x+y;a[j+s]=x-y;}for(int x:a)if(x!=8&&x!=-8)return false;return true;}
int main(int argc,char**argv){
 if(argc!=4)return 2;U begin=std::stoull(argv[1]),count=std::stoull(argv[2]);std::string prefix=argv[3];
 if(begin>(1ULL<<31)||count>(1ULL<<31)-begin)return 3;
 for(int j=0;j<6;++j)for(int x=0;x<64;++x)if(!((x>>j)&1))low[j]|=U(1)<<x;
 Basis w;w.add(~U(0));for(int j=0;j<6;++j)w.add(monomial(1U<<j));std::vector<U> rel;
 for(int j=0;j<8;++j)if(j>=4){w.add(INPUT_REL[j]);rel.push_back(INPUT_REL[j]);}
 if(w.dim!=11)return 4;
 Basis total=w;std::vector<U> qb;
 for(int m=63;m>=0;--m)if(__builtin_popcount(unsigned(m))<=3){U x=monomial(m);if(total.add(x))qb.push_back(x);}
 if(qb.size()!=31||total.dim!=42)return 5;
 U seed=0;for(int x=0;x<64;++x)if((__builtin_popcount(unsigned((x&7)&INPUT_PERM[x>>3]))^((INPUT_H_BITS>>(x>>3))&1))&1)seed|=U(1)<<x;
 if(!bent_derivative(seed)||bent_derivative(0))return 6;
 // Distinct criterion cross-check on 4096 deterministic cubic phases.
 U state=2026111706ULL;
 for(int i=0;i<4096;++i){U f=seed;for(U b:qb){state^=state<<13;state^=state>>7;state^=state<<17;if(state&1)f^=b;}if(bent_derivative(f)!=bent_walsh(f))return 7;}
 std::vector<U> translates{0};for(U r:rel){auto n=translates.size();for(std::size_t j=0;j<n;++j)translates.push_back(translates[j]^r);}
 // Reverse slice order after the zero slice, independent of the primary Gray traversal.
 std::reverse(translates.begin()+1,translates.end());
 U phase=seed,gray=begin^(begin>>1);for(int j=0;j<31;++j)if((gray>>j)&1)phase^=qb[j];
 auto start=std::chrono::steady_clock::now();std::ofstream found(prefix+".found.jsonl");U hits=0;
 for(U i=begin;i<begin+count;++i){bool good=true;for(U t:translates)if(!bent_derivative(phase^t)){good=false;break;}
  if(good){++hits;found<<"{\"index\":"<<i<<",\"relative_phase_hex\":\"0x"<<std::hex<<(phase^seed)<<"\",\"class_normal_form_hex\":\"0x"<<w.reduce(phase^seed)<<std::dec<<"\"}\n";}
  if(i+1<begin+count)phase^=qb[__builtin_ctzll(i+1)];
 }
 double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::ofstream out(prefix+".summary.json");out<<"{\"start\":"<<begin<<",\"count\":"<<count<<",\"hits\":"<<hits<<",\"seconds\":"<<sec<<",\"criterion\":\"all_63_derivatives_balanced\",\"W_rank\":11,\"quotient_rank\":31,\"criterion_crosschecks\":4096}\n";
 std::cout<<"tested="<<count<<" hits="<<hits<<" seconds="<<sec<<"\n";
}
