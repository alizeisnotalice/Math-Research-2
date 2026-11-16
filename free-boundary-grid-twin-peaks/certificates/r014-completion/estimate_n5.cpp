// Lightweight exact state-space estimate for the per-B induced completion DP.
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <algorithm>
using U=__uint128_t;
static std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
int main(){constexpr int n=5,N=25;const uint32_t full=(uint32_t(1)<<N)-1;std::array<uint32_t,N>adj{};for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto d:std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}){int X=x+d[0],Y=y+d[1];if(X>=0&&X<n&&Y>=0&&Y<n)adj[v]|=uint32_t(1)<<(Y*n+X);}}
 auto conn=[&](uint32_t S){if(!S)return false;uint32_t seen=S&-S,fr=seen;while(fr){uint32_t b=fr&-fr;fr-=b;uint32_t q=adj[__builtin_ctz(b)]&S&~seen;seen|=q;fr|=q;}return seen==S;};
 std::array<uint64_t,26>cb{},both{};U allcost=0,bothcost=0;for(uint32_t B=1;B<=full;B++){if(!conn(B))continue;int b=__builtin_popcount(B);U cost=U(1)<<(N-b);cb[b]++;allcost+=cost;if(conn(full^B)){both[b]++;bothcost+=cost;}}
 std::cout<<"{\"n\":5,\"connected_B_cost\":\""<<dec(allcost)<<"\",\"connected_B_and_H_cost\":\""<<dec(bothcost)<<"\",\"by_b\":[";bool f=true;for(int b=1;b<=N;b++)if(cb[b]){if(!f)std::cout<<',';f=false;std::cout<<"{\"b\":"<<b<<",\"connected_B\":"<<cb[b]<<",\"H_connected\":"<<both[b]<<"}";}std::cout<<"]}\n";}
