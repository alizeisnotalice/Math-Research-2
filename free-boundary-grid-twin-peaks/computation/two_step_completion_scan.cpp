// Exact scan of p(S)=sum_{x in dS,y in N(x)\closed(S)} h(Sxy)/h(S)
// on n-by-n free grids, n=4 or 5.  All arithmetic is unsigned __int128.
// h(S) is the connected-completion count h(full)=1,
// h(S)=sum_{x in boundary(S)}h(S union {x}).
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

using Mask=uint64_t;
using U128=unsigned __int128;
constexpr int THREADS=12;
static int n,N;
static Mask FULL;
static std::array<Mask,25> adj;
static U128* h=nullptr;

static int pop(Mask x){return __builtin_popcountll(x);}
static int bit(Mask x){return __builtin_ctzll(x);}
static std::string dec(U128 x){if(!x)return "0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
static Mask neigh(Mask s){Mask r=0;while(s){Mask b=s&-s;s-=b;r|=adj[bit(b)];}return r;}
static bool connected(Mask s){if(!s)return false;Mask seen=s&-s,fr=seen;while(fr){Mask q=fr;fr=0;while(q){Mask b=q&-q;q-=b;fr|=adj[bit(b)];}fr&=s&~seen;seen|=fr;}return seen==s;}

// Exact comparison a/b vs c/d for positive U128 quantities, avoiding the
// potentially 256-bit cross product via continued fractions.
static int cmp_frac(U128 a,U128 b,U128 c,U128 d){
  int sign=1;
  for(;;){U128 qa=a/b,qc=c/d;if(qa!=qc)return qa<qc?-sign:sign;U128 ra=a%b,rc=c%d;
    if(!ra||!rc){if(!ra&&!rc)return 0;return !ra?-sign:sign;}
    a=b;b=ra;c=d;d=rc;sign=-sign;
  }
}
struct Candidate {bool set=false;Mask S=0;int B=0;U128 num=0,den=1;};
static bool better(const Candidate&a,const Candidate&b){return !b.set||(a.set&&cmp_frac(a.num*U128(a.B),a.den,b.num*U128(b.B),b.den)>0);}
static std::string cells(Mask s){std::string z="[";bool first=true;for(int v=0;v<N;v++)if(s&(Mask(1)<<v)){if(!first)z+=",";first=false;z+="["+std::to_string(v%n)+","+std::to_string(v/n)+"]";}return z+"]";}

int main(int argc,char**argv){
 try{
  if(argc!=2)return 2;n=std::atoi(argv[1]);if(n<2||n>5)return 2;N=n*n;FULL=(Mask(1)<<N)-1;
  for(int v=0;v<N;v++){int x=v%n,y=v/n;for(auto [X,Y]:{std::pair{x-1,y},std::pair{x+1,y},std::pair{x,y-1},std::pair{x,y+1}})if(0<=X&&X<n&&0<=Y&&Y<n)adj[v]|=Mask(1)<<(Y*n+X);}
  size_t sz=size_t(FULL)+1;h=(U128*)std::calloc(sz,sizeof(U128));if(!h)throw std::bad_alloc();h[FULL]=1;
  // Dependencies are strictly in the next cardinality layer.  Each layer is
  // partitioned into 12 contiguous numeric ranges, so writes are disjoint.
  for(int k=N-1;k>=1;k--){std::vector<std::thread> th;for(int t=0;t<THREADS;t++)th.emplace_back([&,t]{Mask lo=(FULL*Mask(t))/THREADS+1,hi=(FULL*Mask(t+1))/THREADS+1;for(Mask S=lo;S<hi;S++){if(pop(S)!=k)continue;Mask bd=neigh(S)&~S;U128 z=0;while(bd){Mask b=bd&-bd;bd-=b;z+=h[S|b];}h[S]=z;}});for(auto&t:th)t.join();}
  struct Local {uint64_t tested=0,violations=0;Candidate best,first;};std::array<Local,THREADS> local{};std::vector<std::thread> th;
  for(int t=0;t<THREADS;t++)th.emplace_back([&,t]{Mask lo=(FULL*Mask(t))/THREADS+1,hi=(FULL*Mask(t+1))/THREADS+1;auto&out=local[t];for(Mask S=lo;S<hi;S++){if(!connected(S))continue;Mask bd=neigh(S)&~S;int B=pop(bd);if(!B)continue;out.tested++;Mask closed=S|neigh(S);U128 num=0;for(Mask q=bd;q;){Mask xb=q&-q;q-=xb;Mask ys=adj[bit(xb)]&~closed;while(ys){Mask yb=ys&-ys;ys-=yb;num+=h[S|xb|yb];}}
      Candidate c{true,S,B,num,h[S]};if(better(c,out.best))out.best=c;if(num*U128(B)>U128(4)*h[S]){out.violations++;if(!out.first.set)out.first=c;}
  }});for(auto&t:th)t.join();
  uint64_t tested=0,viol=0;Candidate best,first;for(auto&o:local){tested+=o.tested;viol+=o.violations;if(better(o.best,best))best=o.best;if(!first.set&&o.first.set)first=o.first;}
  U128 A1=0;for(int v=0;v<N;v++)A1+=h[Mask(1)<<v];
  std::cout<<"{\"n\":"<<n<<",\"threads\":"<<THREADS<<",\"A1\":\""<<dec(A1)<<"\",\"tested_connected_prefixes\":"<<tested<<",\"candidate_p_le_4_over_boundary_violations\":"<<viol<<",\"max_p_times_boundary\":{\"numerator\":\""<<dec(best.num*U128(best.B))<<"\",\"denominator\":\""<<dec(best.den)<<"\",\"boundary\":"<<best.B<<",\"S_cells\":"<<cells(best.S)<<"}";
  if(first.set)std::cout<<",\"first_violation\":{\"numerator\":\""<<dec(first.num)<<"\",\"denominator\":\""<<dec(first.den)<<"\",\"boundary\":"<<first.B<<",\"S_cells\":"<<cells(first.S)<<"}";else std::cout<<",\"first_violation\":null";
  std::cout<<"}\n";std::free(h);return 0;
 }catch(const std::exception&e){std::cerr<<"two_step_completion_scan: "<<e.what()<<"\n";std::free(h);return 1;}
}
