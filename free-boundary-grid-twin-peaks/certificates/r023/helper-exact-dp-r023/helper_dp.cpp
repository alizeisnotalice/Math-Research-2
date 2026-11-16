// Exact subset-DP for g/h on K=3 peak sets with a 2+1 domino cover.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
using U=__uint128_t; using M=uint32_t;
std::string ss(U x){if(!x)return"0";std::string s;while(x){s.push_back('0'+x%10);x/=10;}std::reverse(s.begin(),s.end());return s;}
struct R{U words=0;std::vector<U> num;};
int n,N; std::vector<M> nb;
bool ind(M p){for(int v=0;v<N;v++)if((p>>v&1)&&(nb[v]&p))return false;return true;}
R dp(M P,int u,int v){
  M F=(M(1)<<N)-1; std::vector<U> H(F+1),C(F+1),Uu(F+1),Uv(F+1); H[F]=1;
  auto ok=[&](M s,int w){bool e=nb[w]&s;return ((P>>w&1)?!e:e);};
  for(M s=F;s--;){for(int w=0;w<N;w++)if(!(s>>w&1)&&ok(s,w))H[s]+=H[s|(M(1)<<w)];}
  C[0]=1;R z;z.num.assign(N,U(0));
  for(M s=0;s<F;s++)if(C[s]){int k=__builtin_popcount(s);for(int w=0;w<N;w++)if(!(s>>w&1)&&ok(s,w)){
    M t=s|(M(1)<<w);C[t]+=C[s];Uu[t]+=Uu[s]+(w==u?U(k)*C[s]:0);Uv[t]+=Uv[s]+(w==v?U(k)*C[s]:0);
    if(w==v&&(s>>u&1)&&k<N-1)z.num[k]+=(U(k)*C[s]-Uu[s])*H[t];
    if(w==u&&(s>>v&1)&&k<N-1)z.num[k]+=(U(k)*C[s]-Uv[s])*H[t];
  }}
  z.words=C[F];return z;
}
struct Cov{M P;int a,b,u,v,q;};
std::vector<Cov> covers(){
 std::vector<Cov> z;std::set<std::pair<M,std::pair<int,int>>> seen;
 for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)if(nb[a]>>b&1){
  M ca=nb[a]&~(M(1)<<b),cb=nb[b]&~(M(1)<<a),h=ca|cb;
  for(M p=0;p<(M(1)<<N);p++)if(__builtin_popcount(p)==3&&(p&~h)==0&&ind(p)){
   M common=(M(1)<<N)-1;for(int x=0;x<N;x++)if(p>>x&1)common&=nb[x];if(common)continue;
   auto addcover=[&](int h,int other,M heavy,M light){
    int x[3],k=0;for(int w=0;w<N;w++)if(p>>w&1)x[k++]=w;int u=-1,v=-1,q=-1;
    for(int w:x)if(heavy>>w&1){if(u<0)u=w;else v=w;}else q=w;
    assert(u>=0&&v>=0&&q>=0);z.push_back({p,h,other,u,v,q});
   };
   if(__builtin_popcount(p&ca)==2&&__builtin_popcount(p&cb)==1)addcover(a,b,ca,cb);
   if(__builtin_popcount(p&cb)==2&&__builtin_popcount(p&ca)==1)addcover(b,a,cb,ca);
  }
 }return z;
}
int main(int ac,char**av){if(ac!=2&&ac!=4&&ac!=5)return 2;n=std::atoi(av[1]);if(n<3||n>5)return 2;N=n*n;assert(N<32);nb.assign(N,0);
 for(int y=0;y<n;y++)for(int x=0;x<n;x++){int w=y*n+x;for(auto d:{std::pair{-1,0},{1,0},{0,-1},{0,1}}){int X=x+d.first,Y=y+d.second;if(X>=0&&X<n&&Y>=0&&Y<n)nb[w]|=M(1)<<(Y*n+X);}}
 auto cs=covers();if(ac>=4){int a=std::atoi(av[2]),b=std::atoi(av[3]);std::vector<Cov> q;for(auto c:cs)if(c.a==a&&c.b==b&&(ac==4||c.P==M(std::strtoul(av[4],nullptr,10))))q.push_back(c);if(q.empty())return 3;q.resize(1);cs.swap(q);}std::map<std::pair<M,int>,R> cache;U total=0;std::vector<U> num(N);std::set<M> ps;
 for(auto c:cs){int paircode=std::min(c.u,c.v)*N+std::max(c.u,c.v);auto key=std::make_pair(c.P,paircode);if(!cache.count(key))cache[key]=dp(c.P,c.u,c.v);auto&r=cache[key];ps.insert(c.P);total+=r.words;for(int k=0;k<N;k++)num[k]+=r.num[k];}
 std::cout<<"{\"n\":"<<n<<",\"cover_records\":"<<cs.size()<<",\"distinct_P\":"<<ps.size()<<",\"pair_DPs\":"<<cache.size()<<",\"cover_word_total\":\""<<ss(total)<<"\",\"rank_numerators\":{";
 bool f=1;for(int k=1;k<N-1;k++)if(num[k]){if(!f)std::cout<<',';f=0;std::cout<<'\"'<<k<<"\":\""<<ss(num[k])<<"\"";}std::cout<<"},\"records\":[";
 for(size_t i=0;i<cs.size();i++){auto c=cs[i];int pc=std::min(c.u,c.v)*N+std::max(c.u,c.v);auto&r=cache[{c.P,pc}];if(i)std::cout<<',';std::cout<<"{\"P\":"<<c.P<<",\"edge\":["<<c.a<<','<<c.b<<"],\"helper\":"<<c.a<<",\"same\":["<<c.u<<','<<c.v<<"],\"q\":"<<c.q<<",\"words\":\""<<ss(r.words)<<"\"}";}std::cout<<"]}\n";
}
