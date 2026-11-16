// Exact first-merge decomposition with one-peak root moments on P_n square.
// n <= 5.  All stored and accumulated quantities are nonnegative __uint128_t.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using U = __uint128_t;

static std::string dec(U x) {
  if (!x) return "0";
  std::string s;
  while (x) { s.push_back(char('0' + x % 10)); x /= 10; }
  std::reverse(s.begin(), s.end());
  return s;
}

static uint64_t choose64(int n, int k) {
  if (k > n-k) k=n-k;
  uint64_t z=1;
  for (int i=1;i<=k;i++) z=z*(n-k+i)/i;
  return z;
}

static uint32_t nextComb(uint32_t x) {
  uint32_t c=x&-x, r=x+c;
  return (((r^x)>>2)/c)|r;
}

static std::vector<uint32_t> masksOfSize(int n, int k) {
  const int N=n*n;
  std::vector<uint32_t> z;
  z.reserve(size_t(choose64(N,k)));
  if (k==N) { z.push_back((uint32_t(1)<<N)-1); return z; }
  uint32_t x=(uint32_t(1)<<k)-1;
  const uint32_t cap=uint32_t(1)<<N;
  while (x<cap) { z.push_back(x); x=nextComb(x); }
  return z;
}

template<class F>
static void parallelRanges(size_t count, int workers, F f) {
  int used=std::min<int>(workers, int(count));
  std::vector<std::thread> ts;
  ts.reserve(used);
  for (int t=0;t<used;t++) {
    size_t lo=count*size_t(t)/used, hi=count*size_t(t+1)/used;
    ts.emplace_back([=,&f](){ f(lo,hi,t); });
  }
  for (auto& t:ts) t.join();
}

struct Agg { U w=0,d=0; };
struct Top {
  U w=0,d=0, aA=0,aB=0, mxA=0,myA=0,m2A=0,mxB=0,myB=0,m2B=0,h=0;
  uint32_t S=0,A=0,B=0;
  uint8_t s=0, branch=0, span=0, wvertex=0;
};
static void keepTop(std::vector<Top>& z, const Top& x) {
  constexpr size_t LIMIT=20;
  if (z.size()<LIMIT) { z.push_back(x); }
  else {
    auto it=std::min_element(z.begin(),z.end(),[](const Top&a,const Top&b){return a.w<b.w;});
    if (it->w>=x.w) return;
    *it=x;
  }
}

static int branchSpan(uint32_t A, uint32_t B, int n) {
  // max_{u in A,v in B} dist_1(u,v), via the four linear functionals.
  int minA[4]={100,100,100,100}, maxA[4]={-100,-100,-100,-100};
  int minB[4]={100,100,100,100}, maxB[4]={-100,-100,-100,-100};
  const int sx[4]={1,1,-1,-1}, sy[4]={1,-1,1,-1};
  auto scan=[&](uint32_t m,int*lo,int*hi){
    while(m) { uint32_t b=m&-m; m-=b; int v=__builtin_ctz(b),x=v%n,y=v/n;
      for(int j=0;j<4;j++){int q=sx[j]*x+sy[j]*y;lo[j]=std::min(lo[j],q);hi[j]=std::max(hi[j],q);}
    }
  };
  scan(A,minA,maxA); scan(B,minB,maxB);
  int ans=0;
  for(int j=0;j<4;j++) ans=std::max(ans,maxA[j]-minB[j]);
  return ans;
}

static uint32_t component(uint32_t seed, uint32_t allowed,
                          const std::array<uint32_t,25>& adj) {
  uint32_t seen=seed, front=seed;
  while(front) { uint32_t b=front&-front; front-=b; int v=__builtin_ctz(b);
    uint32_t add=adj[v]&allowed&~seen; seen|=add; front|=add;
  }
  return seen;
}

int main(int argc,char**argv) {
  if (argc!=4) return 2;
  const int n=std::stoi(argv[1]), N=n*n, workers=std::stoi(argv[3]);
  if (n<2 || n>5 || workers<1) return 2;
  const uint32_t cap=uint32_t(1)<<N, full=cap-1;
  std::array<uint32_t,25> adj{};
  for(int y=0;y<n;y++) for(int x=0;x<n;x++) {
    int v=y*n+x;
    for (auto d: std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}) {
      int X=x+d[0],Y=y+d[1]; if(X>=0&&X<n&&Y>=0&&Y<n) adj[v]|=uint32_t(1)<<(Y*n+X);
    }
  }
  auto alloc=[&](){ return static_cast<U*>(std::calloc(cap,sizeof(U))); };
  U *a=alloc(), *mx=alloc(), *my=alloc(), *m2=alloc(), *h=alloc();
  if(!a||!mx||!my||!m2||!h) { std::cerr<<"allocation failure\n"; return 3; }
  const auto began=std::chrono::steady_clock::now();
  for(int v=0;v<N;v++) { uint32_t S=uint32_t(1)<<v; int x=v%n,y=v/n;
    a[S]=1; mx[S]=x; my[S]=y; m2[S]=U(x*x+y*y);
  }

  // a and its three root-coordinate moment sums, rank by rank.
  for(int k=2;k<=N;k++) {
    auto masks=masksOfSize(n,k);
    parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){
      for(size_t i=lo;i<hi;i++) { uint32_t S=masks[i]; U za=0,zx=0,zy=0,z2=0;
        for(uint32_t q=S;q;q&=q-1) { uint32_t bit=q&-q; int v=__builtin_ctz(bit); uint32_t T=S^bit;
          if(adj[v]&T) { za+=a[T]; zx+=mx[T]; zy+=my[T]; z2+=m2[T]; }
        }
        a[S]=za; mx[S]=zx; my[S]=zy; m2[S]=z2;
      }
    });
  }

  // h(S) is a no-new-birth completion count.  Supersets have one larger rank.
  h[full]=1;
  for(int k=N-1;k>=1;k--) {
    auto masks=masksOfSize(n,k);
    parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){
      for(size_t i=lo;i<hi;i++) { uint32_t S=masks[i]; U z=0, out=full^S;
        for(uint32_t q=out;q;q&=q-1) { uint32_t bit=q&-q; int v=__builtin_ctz(bit); if(adj[v]&S) z+=h[S|bit]; }
        h[S]=z;
      }
    });
  }

  const int MAXMIN=N/2, MAXSPAN=2*(n-1);
  const size_t CELLS=size_t(N+1)*(MAXMIN+1)*(MAXSPAN+1);
  std::vector<std::vector<Agg>> local{size_t(workers),std::vector<Agg>(CELLS)};
  std::vector<std::vector<Top>> localTop{size_t(workers)};
  std::vector<uint64_t> localStates(size_t(workers),0),localEvents(size_t(workers),0);
  auto cell=[&](int s,int mn,int span){ return (size_t(s)*(MAXMIN+1)+mn)*(MAXSPAN+1)+span; };

  for(int s=3;s<=N;s++) {
    auto masks=masksOfSize(n,s);
    parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int tid){
      auto& ledger=local[tid]; auto& tops=localTop[tid];
      uint64_t states=0,events=0;
      for(size_t i=lo;i<hi;i++) { uint32_t S=masks[i];
        if(!h[S] || !a[S]) continue; // a(S)>0 iff G[S] is connected.
        ++states;
        for(uint32_t q=S;q;q&=q-1) { uint32_t bit=q&-q; int w=__builtin_ctz(bit); uint32_t T=S^bit;
          if(!T) continue;
          uint32_t A=component(T&-T,T,adj), B=T^A;
          if(!B || component(B&-B,B,adj)!=B) continue;
          int as=__builtin_popcount(A), bs=__builtin_popcount(B), mn=std::min(as,bs);
          int span=branchSpan(A,B,n);
          U comb=choose64(s-1,as);
          U ww=h[S]*comb*a[A]*a[B];
          U left=a[B]*m2[A]+a[A]*m2[B];
          U right=U(2)*(mx[A]*mx[B]+my[A]*my[B]);
          if(left<right) { std::cerr<<"negative squared-distance numerator\n"; std::abort(); }
          U dd=h[S]*comb*(left-right);
          Agg& g=ledger[cell(s,mn,span)]; g.w+=ww; g.d+=dd; ++events;
          Top rec; rec.w=ww;rec.d=dd;rec.aA=a[A];rec.aB=a[B];rec.mxA=mx[A];rec.myA=my[A];rec.m2A=m2[A];
          rec.mxB=mx[B];rec.myB=my[B];rec.m2B=m2[B];rec.h=h[S];rec.S=S;rec.A=A;rec.B=B;
          rec.s=s;rec.branch=mn;rec.span=span;rec.wvertex=w; keepTop(tops,rec);
        }
      }
      localStates[tid]+=states; localEvents[tid]+=events;
    });
  }

  std::vector<Agg> ledger(CELLS); U totalW=0,totalD=0;
  uint64_t states=0,events=0;
  std::vector<Top> tops;
  for(int t=0;t<workers;t++) { states+=localStates[t];events+=localEvents[t];
    for(size_t i=0;i<CELLS;i++){ledger[i].w+=local[t][i].w;ledger[i].d+=local[t][i].d;}
    for(const auto&r:localTop[t]) keepTop(tops,r);
  }
  for(const auto&g:ledger){totalW+=g.w;totalD+=g.d;}
  std::sort(tops.begin(),tops.end(),[](const Top&x,const Top&y){return x.w>y.w;});

  std::ofstream o(argv[2]);
  o<<"{\n  \"n\": "<<n<<",\n  \"workers\": "<<workers
   <<",\n  \"state_array_bytes\": "<<uint64_t(5)*cap*sizeof(U)
   <<",\n  \"connected_S_examined\": "<<states<<",\n  \"first_merge_events\": "<<events
   <<",\n  \"A2_from_W\": \""<<dec(totalW)<<"\",\n  \"D_euclidean_squared_from_W\": \""<<dec(totalD)<<"\",\n"
   <<"  \"rank_minbranch_span\": [\n";
  bool first=true;
  for(int s=3;s<=N;s++) for(int mn=1;mn<=MAXMIN;mn++) for(int sp=0;sp<=MAXSPAN;sp++) {
    const Agg&g=ledger[cell(s,mn,sp)]; if(!g.w) continue;
    if(!first)o<<",\n"; first=false;
    o<<"    {\"s\": "<<s<<", \"min_branch\": "<<mn<<", \"branch_span_manhattan\": "<<sp
     <<", \"W\": \""<<dec(g.w)<<"\", \"D_euclidean_squared\": \""<<dec(g.d)<<"\"}";
  }
  o<<"\n  ],\n  \"top_W_states\": [\n";
  for(size_t i=0;i<tops.size();i++) { const Top&r=tops[i]; if(i)o<<",\n";
    o<<"    {\"S_mask\": "<<r.S<<", \"w\": "<<int(r.wvertex)<<", \"A_mask\": "<<r.A<<", \"B_mask\": "<<r.B
     <<", \"s\": "<<int(r.s)<<", \"min_branch\": "<<int(r.branch)<<", \"branch_span_manhattan\": "<<int(r.span)
     <<", \"h\": \""<<dec(r.h)<<"\", \"aA\": \""<<dec(r.aA)<<"\", \"aB\": \""<<dec(r.aB)
     <<"\", \"MxA\": \""<<dec(r.mxA)<<"\", \"MyA\": \""<<dec(r.myA)<<"\", \"M2A\": \""<<dec(r.m2A)
     <<"\", \"MxB\": \""<<dec(r.mxB)<<"\", \"MyB\": \""<<dec(r.myB)<<"\", \"M2B\": \""<<dec(r.m2B)
     <<"\", \"W\": \""<<dec(r.w)<<"\", \"D_euclidean_squared\": \""<<dec(r.d)<<"\"}";
  }
  o<<"\n  ],\n  \"wall_seconds\": "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<"\n}\n";
  o.close();
  std::cout<<"n="<<n<<" A2="<<dec(totalW)<<" D="<<dec(totalD)<<" states="<<states<<" events="<<events<<"\n";
  std::free(a);std::free(mx);std::free(my);std::free(m2);std::free(h);
}
