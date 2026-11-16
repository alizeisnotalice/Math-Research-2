#include <cstdint>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <cassert>
using U=__uint128_t;
std::string dec(U z){std::string s;do{s.push_back('0'+z%10);z/=10;}while(z);std::reverse(s.begin(),s.end());return s;}
int main(int argc,char**argv){assert(argc==2);int n=std::stoi(argv[1]);assert(n>=3&&n<=5);int N=n*n;uint32_t all=(1u<<N)-1,l=0,r=0;for(int y=0;y<n;y++){l|=1u<<(y*n);r|=1u<<(y*n+n-1);}std::vector<U> h(size_t(all)+1);h[all]=1;
for(uint32_t s=all;s--;){uint32_t q=(((s<<1)&~l)|((s>>1)&~r)|(s<<n)|(s>>n))&all&~s;while(q){uint32_t v=q&-q;q-=v;h[s]+=h[s|v];}}
std::cout<<"{\"n\":"<<n<<",\"counts\":{\"1\":\""<<dec(h[1])<<"\",\""<<(1u<<(n+1))<<"\":\""<<dec(h[1u<<(n+1)])<<"\"";
for(int a=1;a<=std::min(3,n-2);a++){uint32_t m=1;for(int y=1;y<=a;y++)for(int x=1;x<=a;x++)m|=1u<<(y*n+x);std::cout<<",\""<<m<<"\":\""<<dec(h[m])<<"\",\""<<(m|2)<<"\":\""<<dec(h[m|2])<<"\"";}
std::cout<<"}}\n";}
