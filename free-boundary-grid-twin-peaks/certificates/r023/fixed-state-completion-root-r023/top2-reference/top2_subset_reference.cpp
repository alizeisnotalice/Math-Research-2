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

U B=0,D=0,Fm=0;__int128 Gm=0;std::vector<U> hist(2*n-1);
for(int u=0;u<N;u++)for(int v=u+1;v<N;v++){
 int dx=u%n-v%n,dy=u/n-v/n,d=std::abs(dx)+std::abs(dy);if(d<2)continue;
 U w=2*h[(1u<<u)|(1u<<v)];B+=w;hist[d]+=w;D+=w*(dx*dx+dy*dy);
 int ux=2*(u%n)-(n-1),uy=2*(u/n)-(n-1),vx=2*(v%n)-(n-1),vy=2*(v/n)-(n-1);
 Fm+=w*(ux*ux+uy*uy+vx*vx+vy*vy);Gm+=(__int128)w*(ux*vx+uy*vy);
}
assert((__int128)Fm-2*Gm==4*(__int128)D);
std::cout<<"{\"n\":"<<n<<",\"Btop2\":\""<<dec(B)<<"\",\"D\":\""<<dec(D)<<"\",\"F\":\""<<dec(Fm)<<"\",\"G\":\""<<(Gm<0?"-":"")<<dec(Gm<0?-Gm:Gm)<<"\",\"hist\":[";
for(int d=0;d<2*n-1;d++){if(d)std::cout<<',';std::cout<<'\"'<<dec(hist[d])<<'\"';}std::cout<<"]}\n";
}
