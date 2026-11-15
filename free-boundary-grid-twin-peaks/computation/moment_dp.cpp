// Exact birth-count and peak squared-distance DP. N<=25 keeps all values <2^128.
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using U = unsigned __int128;
std::string dec(U x) { if (!x) return "0"; std::string s; while(x){s.push_back('0'+x%10);x/=10;}return std::string(s.rbegin(),s.rend());}
struct Cell { U a, sx, sy, ss, b, dd; };
int main(int argc,char**argv){
 if(argc!=2)return 2;
 int n=std::stoi(argv[1]),N=n*n;
 if(n<1||n>5)return 2;
 uint32_t size=1u<<N;
 auto start=std::chrono::steady_clock::now();
 Cell* f=(Cell*)std::calloc(size,sizeof(Cell));
 if(!f)return 3;
 std::vector<uint32_t>adj(N);
 for(int v=0;v<N;v++) for(int u=0;u<N;u++)
  if(std::abs(v%n-u%n)+std::abs(v/n-u/n)==1)adj[v]|=1u<<u;
 for(uint32_t S=1;S<size;S++){
  Cell z{};
  if(!(S&(S-1))){int v=__builtin_ctz(S);z.a=1;z.sx=v%n;z.sy=v/n;z.ss=(v%n)*(v%n)+(v/n)*(v/n);}
  else for(uint32_t bits=S;bits;bits&=bits-1){
   int v=__builtin_ctz(bits);uint32_t T=S^(1u<<v);const Cell&t=f[T];
   if(adj[v]&T){z.a+=t.a;z.sx+=t.sx;z.sy+=t.sy;z.ss+=t.ss;z.b+=t.b;z.dd+=t.dd;}
   else {U x=v%n,y=v/n;z.b+=t.a;z.dd+=t.ss+(x*x+y*y)*t.a-2*(x*t.sx+y*t.sy);}
  }
  f[S]=z;
 }
 const Cell&z=f[size-1];
 std::cout<<"{\"n\":"<<n<<",\"A1\":"<<dec(z.a)<<",\"A2\":"<<dec(z.b)
 <<",\"squared_euclidean_sum\":"<<dec(z.dd)<<",\"Q_numerator\":"<<dec(z.dd)
 <<",\"Q_denominator\":"<<dec(U(N)*z.b)<<",\"allocated_bytes\":"<<uint64_t(size)*sizeof(Cell)
 <<",\"wall_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 std::free(f);
}
