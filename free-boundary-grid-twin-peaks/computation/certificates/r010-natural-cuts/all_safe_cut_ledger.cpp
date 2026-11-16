// Exact ledger for all connected left/right cuts of P_n square P_n.
// A root in an outer column is admitted iff *all* graph neighbours stay in
// its side.  This protects the root from cross-edge suppression.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <gmpxx.h>
using U128=__uint128_t; using Big=mpz_class;
static std::string dec(U128 x){if(!x)return"0";std::string s;while(x){s.push_back('0'+x%10);x/=10;}std::reverse(s.begin(),s.end());return s;}
struct Cut {uint32_t A,Br;uint8_t rootsA,rootsB,size,cross;U128 a=0,b=0;};
struct Work {
 int n,N;uint32_t full;std::array<uint32_t,25> adj{};std::array<std::unordered_map<uint32_t,U128>,5> memo;
 Work(int nn):n(nn),N(nn*nn),full((1u<<N)-1){assert(n>=2&&n<=5);for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto d:{std::pair<int,int>{1,0},{-1,0},{0,1},{0,-1}}){int X=x+d.first,Y=y+d.second;if(X>=0&&X<n&&Y>=0&&Y<n)adj[v]|=1u<<(Y*n+X);}}}
 bool conn(uint32_t s)const{if(!s)return false;uint32_t seen=s&-s,fr=seen;while(fr){uint32_t nb=0,x=fr;while(x){uint32_t b=x&-x;x-=b;nb|=adj[__builtin_ctz(b)];}fr=nb&s&~seen;seen|=fr;}return seen==s;}
 uint32_t refl(uint32_t s)const{uint32_t z=0;for(int y=0;y<n;y++)for(int x=0;x<n;x++)if(s&(1u<<(y*n+x)))z|=1u<<(y*n+n-1-x);return z;}
 uint8_t safeRoots(uint32_t s)const{uint8_t z=0;for(int y=0;y<n;y++){int r=y*n;if((s&(1u<<r)) && (adj[r]&~s)==0)z|=uint8_t(1u<<y);}return z;}
 U128 rootCount(uint32_t s,int y){int r=y*n;if(!(s&(1u<<r)))return 0;auto&M=memo[y];auto it=M.find(s);if(it!=M.end())return it->second;U128 z=0;if(s==(1u<<r))z=1;else{uint32_t x=s&~(1u<<r);while(x){uint32_t b=x&-x;x-=b;if(conn(s^b))z+=rootCount(s^b,y);}}M.emplace(s,z);return z;}
 int cross(uint32_t A)const{uint32_t B=full^A;int z=0;for(int v=0;v<N;v++)if(A&(1u<<v))z+=__builtin_popcount(adj[v]&B);return z;}
 static Big choose(int N,int k){Big z=1;for(int i=1;i<=k;i++)z=z*(N-k+i)/i;return z;}
 static Big big(U128 x){Big z=static_cast<unsigned long>(x>>64);z<<=64;z+=static_cast<unsigned long>(x);return z;}
 void run(const char*out){
   std::vector<Cut> cuts;cuts.reserve(n==5?2500000:10000);
   for(uint32_t A=1;A<=full;A++){uint32_t B=full^A;if(!conn(A)||!conn(B))continue;uint32_t Br=refl(B);uint8_t ra=safeRoots(A),rb=safeRoots(Br);if(!ra||!rb)continue;cuts.push_back({A,Br,ra,rb,uint8_t(__builtin_popcount(A)),uint8_t(cross(A))});}
   for(int y=0;y<n;y++){
      memo[y].reserve(n==5?2000000:100000);
      uint8_t bit=uint8_t(1u<<y);for(auto&c:cuts){if(c.rootsA&bit)c.a+=rootCount(c.A,y);if(c.rootsB&bit)c.b+=rootCount(c.Br,y);}std::cerr<<"root="<<y<<" memo="<<memo[y].size()<<"\n";
   }
   std::map<std::pair<int,int>,Big> by;Big raw=0;int maxe=0,usable=0;for(auto&c:cuts){if(!c.a||!c.b)continue;Big z=choose(N,c.size);z*=big(c.a);z*=big(c.b);by[{c.size,c.cross}]+=z;raw+=z;maxe=std::max(maxe,int(c.cross));usable++;}
   Big w=0;for(auto&[k,z]:by)w+=z<<(maxe-k.second);
   std::ofstream o(out);o<<"{\n  \"n\": "<<n<<",\n  \"connected_safe_cuts\": "<<cuts.size()<<",\n  \"usable_cuts\": "<<usable<<",\n  \"raw_sum\": \""<<raw<<"\",\n  \"weighted_numerator\": \""<<w<<"\",\n  \"common_denominator\": "<<(Big(1)<<maxe)<<",\n  \"max_cross_edges\": "<<maxe<<",\n  \"by_size_cross\": [\n";bool first=true;for(auto&[k,z]:by){if(!first)o<<",\n";first=false;o<<"    {\"A_size\": "<<k.first<<", \"cross_edges\": "<<k.second<<", \"raw_count\": \""<<z<<"\"}";}o<<"\n  ],\n  \"root_memo_states\": [";for(int y=0;y<n;y++){if(y)o<<", ";o<<memo[y].size();}o<<"]\n}\n";
   std::ofstream q(std::string(out).substr(0,std::string(out).rfind('.'))+"-records.json");q<<"{\n  \"n\": "<<n<<",\n  \"records\": [\n";bool rf=true;for(auto&c:cuts){if(!c.a||!c.b)continue;if(!rf)q<<",\n";rf=false;Big z=choose(N,c.size);z*=big(c.a);z*=big(c.b);q<<"    {\"A_mask\": "<<c.A<<", \"A_size\": "<<int(c.size)<<", \"cross_edges\": "<<int(c.cross)<<", \"a_left\": \""<<dec(c.a)<<"\", \"a_right\": \""<<dec(c.b)<<"\", \"raw_count\": \""<<z<<"\"}";}q<<"\n  ]\n}\n";
   std::cerr<<"n="<<n<<" connected_safe="<<cuts.size()<<" usable="<<usable<<" raw="<<raw<<" weighted="<<w<<"/"<<(Big(1)<<maxe)<<"\n";
 }
};
int main(int argc,char**argv){if(argc!=3)return 2;Work(std::stoi(argv[1])).run(argv[2]);}
