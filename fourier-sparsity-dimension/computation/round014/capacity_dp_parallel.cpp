#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>
using U=uint64_t;
template<class F> void parallel(const std::vector<int>& jobs,int nt,F fn){
 std::atomic<size_t> next{0}; std::vector<std::thread> ts;
 for(int t=0;t<nt;t++)ts.emplace_back([&](){for(;;){size_t j=next.fetch_add(16);if(j>=jobs.size())break;for(size_t i=j;i<std::min(j+16,jobs.size());i++)fn(jobs[i]);}});
 for(auto& t:ts)t.join();
}
int main(int argc,char** argv){
 if(argc<4)return 2;int nt=std::stoi(argv[1]),d=std::stoi(argv[2]);
 if(nt<1||nt>64||d<1||d>20||argc!=d+3)return 2;
 auto began=std::chrono::steady_clock::now();int M=1<<d;std::vector<U> cw(M);
 for(int i=0;i<d;i++){U g=std::stoull(argv[i+3],nullptr,0);for(int s=0;s<(1<<i);s++)cw[s|(1<<i)]=cw[s]^g;}
 auto sorted=cw;std::sort(sorted.begin(),sorted.end());if(std::unique(sorted.begin(),sorted.end())!=sorted.end())return 3;
 sorted.clear();sorted.shrink_to_fit();int root=std::find(cw.begin(),cw.end(),~U(0))-cw.begin();if(root==M)return 4;
 std::vector<int> sub(M),dp(M,-1),pick(M),jobs(M);std::iota(jobs.begin(),jobs.end(),0);
 parallel(jobs,nt,[&](int s){int count=0;for(int a=0;a<M;a++)count+=((cw[a]&~cw[s])==0);sub[s]=__builtin_ctz((unsigned)count);});
 std::cerr<<"subcode counts completed\n";
 std::vector<std::vector<int>> levels(65);for(int s=1;s<M;s++)levels[__builtin_popcountll(cw[s])].push_back(s);dp[0]=0;
 for(int w=1;w<=64;w++){
  parallel(levels[w],nt,[&](int s){U R=cw[s],anchor=R&(-R);int best=-1,choice=0;
   for(int a=1;a<M;a++)if((cw[a]&~R)==0 && (cw[a]&anchor)){int v=sub[a]+1+dp[s^a];if(v>best){best=v;choice=a;}}
   dp[s]=best;pick[s]=choice;
  });
  if(!levels[w].empty())std::cerr<<"weight "<<w<<" states "<<levels[w].size()<<" done\n";
 }
 std::cout<<"{\"threads\":"<<nt<<",\"D\":"<<d<<",\"kappa\":"<<dp[root]<<",\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<",\"partition\":[";
 bool first=true;for(int s=root;s;){int a=pick[s];if(!a)return 5;if(!first)std::cout<<",";first=false;std::cout<<"\"0x"<<std::hex<<cw[a]<<std::dec<<"\"";s^=a;}
 std::cout<<"]}\n";
}
