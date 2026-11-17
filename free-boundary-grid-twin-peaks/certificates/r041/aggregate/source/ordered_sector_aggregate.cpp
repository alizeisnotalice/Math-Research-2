// Exact ordered first-merge aggregate DP on P_n, intended for n<=5.
// n=5 allocation is deliberately opt-in: --allow-big-memory is required.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <new>
#include <string>
#include <vector>

using U = unsigned __int128;
static std::string dec(U x) { if (!x) return "0"; std::string s; while(x){s.push_back(char('0'+x%10));x/=10;} std::reverse(s.begin(),s.end());return s; }
static uint64_t choose64(int n,int k) { if(k<0||k>n)return 0;k=std::min(k,n-k);uint64_t z=1;for(int i=1;i<=k;i++)z=z*uint64_t(n-k+i)/uint64_t(i);return z; }
static uint32_t next_comb(uint32_t x) { const uint32_t c=x&-x,r=x+c; return (((r^x)>>2)/c)|r; }
static uint32_t component(uint32_t seed,uint32_t allowed,const std::vector<uint32_t>& adj) {
  uint32_t seen=seed,front=seed; while(front){uint32_t b=front&-front;front^=b;uint32_t add=adj[__builtin_ctz(b)]&allowed&~seen;seen|=add;front|=add;}return seen;
}
struct Cell { U count=0,d2=0; };

int main(int argc,char**argv) {
  if(argc<3){std::cerr<<"usage: ordered_sector_aggregate N OUTPUT [--allow-big-memory] [--plan]\n";return 2;}
  const int n=std::stoi(argv[1]); bool allow=false,plan=false;
  for(int i=3;i<argc;i++){std::string x=argv[i];if(x=="--allow-big-memory")allow=true;else if(x=="--plan")plan=true;else{std::cerr<<"unknown option\n";return 2;}}
  if(n<3||n>5){std::cerr<<"hard guard: only 3 <= n <= 5\n";return 2;}
  const int N=n*n; const uint64_t cap64=uint64_t(1)<<N; const uint32_t cap=uint32_t(cap64),full=cap-1;
  const uint64_t rooted_bytes=cap64*4*sizeof(U), h_bytes=cap64*sizeof(U), unionadj_bytes=cap64*sizeof(uint32_t);
  const uint64_t core_bytes=rooted_bytes+h_bytes+unionadj_bytes;
  std::cout<<"stage=plan n="<<n<<" states="<<cap64<<" rooted_bytes="<<rooted_bytes<<" h_bytes="<<h_bytes<<" unionadj_bytes="<<unionadj_bytes<<" core_bytes="<<core_bytes<<"\n";
  if(plan)return 0;
  if(n==5&&!allow){std::cerr<<"refusing n=5 without --allow-big-memory; run --plan first\n";return 5;}
  if(std::filesystem::exists(argv[2])){std::cerr<<"refusing non-fresh output path\n";return 4;}
  std::vector<uint32_t> adj(N,0);
  for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto d:std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}){int X=x+d[0],Y=y+d[1];if(0<=X&&X<n&&0<=Y&&Y<n)adj[v]|=uint32_t(1)<<(Y*n+X);}}
  try {
    // Four rooted aggregates (a,sx,sy,s2), completion h, and union-neighbor mask.
    std::vector<U> a(cap),sx(cap),sy(cap),s2(cap),h(cap);
    std::vector<uint32_t> uadj(cap,0);
    std::cout<<"stage=allocated core_bytes="<<core_bytes<<"\n";
    for(uint32_t S=1;S<cap;S++){uint32_t bit=S&-S;uadj[S]=uadj[S^bit]|adj[__builtin_ctz(bit)];}
    for(int v=0;v<N;v++){uint32_t S=uint32_t(1)<<v;int x=v%n,y=v/n;a[S]=1;sx[S]=x;sy[S]=y;s2[S]=U(x*x+y*y);}
    for(int k=2;k<=N;k++){
      uint32_t S=(uint32_t(1)<<k)-1; while(S<cap){U za=0,zx=0,zy=0,z2=0;for(uint32_t q=S;q;q&=q-1){uint32_t bit=q&-q,T=S^bit;if(uadj[T]&bit){za+=a[T];zx+=sx[T];zy+=sy[T];z2+=s2[T];}}a[S]=za;sx[S]=zx;sy[S]=zy;s2[S]=z2;
        if(S==full)break;S=next_comb(S);
      }
      if(n==5 && (k==5||k==10||k==15||k==20||k==25))std::cerr<<"stage=rooted_histories rank="<<k<<"\n";
    }
    h[full]=1;
    for(int k=N-1;k>=1;k--){uint32_t S=(uint32_t(1)<<k)-1;while(S<cap){U z=0;for(uint32_t q=(uadj[S]&~S);q;q&=q-1)z+=h[S|(q&-q)];h[S]=z;if(S==full)break;S=next_comb(S);} if(n==5&&(k==20||k==15||k==10||k==5||k==1))std::cerr<<"stage=completions rank="<<k<<"\n";}
    std::vector<Cell> ledger(size_t(N+1)*(N+1)); auto at=[&](int q,int b)->Cell&{return ledger[size_t(q)*(N+1)+b];};
    U A2=0,D2=0; uint64_t states=0,cuts=0;
    std::cerr<<"stage=sector_enumeration\n";
    for(int s=3;s<=N;s++){uint32_t S=(uint32_t(1)<<s)-1;while(S<cap){if(a[S]&&h[S]){++states;for(uint32_t qw=S;qw;qw&=qw-1){uint32_t wbit=qw&-qw;int w=__builtin_ctz(wbit);uint32_t T=S^wbit,B=component(T&-T,T,adj),C=T^B;if(!C||component(C&-C,C,adj)!=C)continue;if(!(adj[w]&B)||!(adj[w]&C))continue;++cuts;const int qrank=s-1,bB=__builtin_popcount(B),bC=__builtin_popcount(C);
          for(int dir=0;dir<2;dir++){uint32_t L=dir?C:B,H=dir?B:C;int b=dir?bC:bB;U inter=choose64(qrank-1,b);U weight=h[S]*inter*a[L]*a[H];U lhs=s2[L]*a[H]+s2[H]*a[L],rhs=U(2)*(sx[L]*sx[H]+sy[L]*sy[H]);if(lhs<rhs){std::cerr<<"negative M_D2\n";return 6;}U d=lhs-rhs;at(qrank,b).count+=weight;at(qrank,b).d2+=h[S]*inter*d;A2+=weight;D2+=h[S]*inter*d;
          }
        }}
        if(S==full)break;S=next_comb(S);
      }
      if(n==5&&(s==5||s==10||s==15||s==20||s==25))std::cerr<<"stage=sector rank="<<s<<"\n";
    }
    std::ofstream out(argv[2]); if(!out){std::cerr<<"cannot create fresh output\n";return 7;}
    out<<"{\n  \"schema\": \"ordered-first-merge-aggregate-v1\",\n  \"n\": "<<n<<",\n  \"workers\": 1,\n  \"manhattan_histogram_included\": false,\n  \"core_allocation_bytes\": "<<core_bytes<<",\n  \"A1\": \""<<dec(a[full])<<"\",\n  \"A2\": \""<<dec(A2)<<"\",\n  \"D_euclidean_squared_sum\": \""<<dec(D2)<<"\",\n  \"eligible_connected_states\": "<<states<<",\n  \"merge_cuts\": "<<cuts<<",\n  \"ordered_q_lowb\": [\n";
    bool first=true;for(int q=2;q<N;q++)for(int b=1;b<q;b++){const Cell&z=at(q,b);if(!z.count)continue;if(!first)out<<",\n";first=false;out<<"    {\"q\": "<<q<<", \"low_b\": "<<b<<", \"count\": \""<<dec(z.count)<<"\", \"D_euclidean_squared_sum\": \""<<dec(z.d2)<<"\"}";}
    out<<"\n  ]\n}\n";std::cout<<"stage=done n="<<n<<" A1="<<dec(a[full])<<" A2="<<dec(A2)<<" D2="<<dec(D2)<<"\n";
  } catch(const std::bad_alloc&){std::cerr<<"allocation failure after plan; no output committed\n";return 8;}
  return 0;
}
