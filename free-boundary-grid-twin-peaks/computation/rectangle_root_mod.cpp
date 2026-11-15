// Exact unique-peak root-band count by two coprime machine-word moduli.
// The two residues reconstruct the integer because it is at most (rc)!.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using U=unsigned __int128;
std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}return {s.rbegin(),s.rend()};}
constexpr uint64_t M=UINT64_MAX; // 2^64-1, coprime to 2^64.
uint64_t add_mod_m(uint64_t a,uint64_t b){return uint64_t((U(a)+b)%M);}
uint64_t run(int rows,int cols,uint64_t allowed){
 int N=rows*cols;uint32_t size=uint32_t(1)<<N,full=size-1;std::vector<uint32_t>adj(N);
 for(int v=0;v<N;v++)for(int w=0;w<N;w++)if(std::abs(v%cols-w%cols)+std::abs(v/cols-w/cols)==1)adj[v]|=uint32_t(1)<<w;
 std::vector<uint64_t> f(size);
 for(uint32_t S=1;S<size;S++){
  if(!(S&(S-1))){int v=__builtin_ctz(S);f[S]=(allowed>>v)&1;continue;}
  uint64_t z=0;for(uint32_t bits=S;bits;bits&=bits-1){int v=__builtin_ctz(bits);uint32_t T=S^(uint32_t(1)<<v);if(adj[v]&T)z=add_mod_m(z,f[T]);}f[S]=z;
 }
 return f[full];
}
int main(int argc,char**argv){
 if(argc!=3)return 2;int rows=std::stoi(argv[1]),cols=std::stoi(argv[2]),N=rows*cols;if(N>24||rows<1||cols<2)return 2;
 uint64_t allowed=0;for(int y=0;y<rows;y++)for(int x=0;x<cols;x++)if(x<=(cols-1)/2)allowed|=uint64_t(1)<<(y*cols+x);
 // First pass modulo 2^64: unsigned overflow is the modulus.
 int NN=N;uint32_t size=uint32_t(1)<<N;std::vector<uint32_t>adj(NN);for(int v=0;v<NN;v++)for(int w=0;w<NN;w++)if(std::abs(v%cols-w%cols)+std::abs(v/cols-w/cols)==1)adj[v]|=uint32_t(1)<<w;
 std::vector<uint64_t> f(size);for(uint32_t S=1;S<size;S++){if(!(S&(S-1))){int v=__builtin_ctz(S);f[S]=(allowed>>v)&1;continue;}uint64_t z=0;for(uint32_t bits=S;bits;bits&=bits-1){int v=__builtin_ctz(bits),T=S^(uint32_t(1)<<v);if(adj[v]&T)z+=f[T];}f[S]=z;}uint64_t r0=f.back();f.clear();f.shrink_to_fit();
 uint64_t r1=run(rows,cols,allowed);
 // x=r0+2^64*t and 2^64=1 mod M, so t=r1-r0 mod M.
 uint64_t t=(r1>=r0)?r1-r0:uint64_t(U(r1)+M-r0);U x=U(r0)+(U(1)<<64)*t;
 std::cout<<"rows="<<rows<<" cols="<<cols<<" L="<<dec(x)<<" residue2_64="<<r0<<" residue2_64_minus_1="<<r1<<"\n";
}
