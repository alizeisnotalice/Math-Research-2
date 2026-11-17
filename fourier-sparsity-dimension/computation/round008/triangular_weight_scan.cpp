// Exhaust all 2^31 triangular masks on F2^8, k=4. Integer computation.
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
using U=std::uint64_t;
using Mask=std::array<U,4>;
int main(int argc,char**argv){
 if(argc!=4)return 2;
 U begin=std::stoull(argv[1]),count=std::stoull(argv[2]);
 if(begin>(1ULL<<31)||count>(1ULL<<31)-begin)return 3;
 std::array<Mask,31> basis{};
 for(unsigned v=0;v<16;++v)for(unsigned u=0;u<16;++u){
   unsigned x=u*16+v;basis[v][x/64]|=1ULL<<(x%64);
 }
 unsigned at=16;
 for(unsigned j=0;j<4;++j)for(unsigned p=0;p<(1U<<j);++p){
   for(unsigned v=0;v<16;++v)for(unsigned u=0;u<16;++u)
     if(((u>>j)&1)&&((v&((1U<<j)-1))==p)){
       unsigned x=u*16+v;basis[at][x/64]|=1ULL<<(x%64);
     }
   ++at;
 }
 if(at!=31)return 4;
 U gray=begin^(begin>>1);Mask word{};
 for(unsigned j=0;j<31;++j)if((gray>>j)&1)
   for(unsigned b=0;b<4;++b)word[b]^=basis[j][b];
 std::array<U,257> histogram{},activity{};
 auto start=std::chrono::steady_clock::now();
 for(U index=begin;index<begin+count;++index){
   unsigned wt=0;for(U b:word)wt+=__builtin_popcountll(b);
   ++histogram[wt];activity[wt]|=gray;
   if(index+1<begin+count){
     unsigned change=__builtin_ctzll(index+1);gray^=1ULL<<change;
     for(unsigned b=0;b<4;++b)word[b]^=basis[change][b];
   }
 }
 double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::ofstream out(argv[3]);
 out<<"{\"status\":\"passed\",\"k\":4,\"q\":16,\"dimV\":31,\"start\":"<<begin
    <<",\"count\":"<<count<<",\"elapsed_seconds\":"<<elapsed<<",\"histogram\":[";
 for(unsigned i=0;i<257;++i){if(i)out<<",";out<<histogram[i];}
 out<<"],\"parameter_activity_by_weight\":[";
 for(unsigned i=0;i<257;++i){if(i)out<<",";out<<activity[i];}
 out<<"]}\n";
 std::cout<<"tested="<<count<<" seconds="<<elapsed<<"\n";
}
