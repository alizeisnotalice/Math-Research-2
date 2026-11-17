// Independent low-memory RM(3,6) check: direct ANF weights and coset collisions.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>
using U=std::uint64_t;
static int pop(U x){return __builtin_popcountll(x);}
static U monomial(std::initializer_list<int> v){U out=0;for(int x=0;x<64;++x){int b=1;for(int i:v)b&=(x>>i)&1;if(b)out|=1ULL<<x;}return out;}
static std::string hx(U x){std::ostringstream o;o<<"0x"<<std::hex<<x;return o.str();}
struct K{int a,b,c;};
static std::vector<K> make_kernels(){std::set<int>s;std::vector<K>o;for(int a=1;a<64;++a)for(int b=a+1;b<64;++b){int c=a^b;std::array<int,3>q{a,b,c};std::sort(q.begin(),q.end());int key=q[0]|q[1]<<6|q[2]<<12;if(s.insert(key).second)o.push_back({q[0],q[1],q[2]});}if(o.size()!=651)throw std::runtime_error("kernel count");return o;}
static bool collision(U support,const K&k,int* witness=nullptr){for(int x=0;x<64;++x)if((support>>x)&1)for(int h:{k.a,k.b,k.c})if((support>>(x^h))&1){if(witness)*witness=h;return true;}return false;}
static int first_injective(U support,const std::vector<K>&ks,std::uint64_t* tests=nullptr){for(int i=0;i<(int)ks.size();++i){if(tests)++*tests;if(!collision(support,ks[i]))return i;}return -1;}
struct Stat{std::atomic<U> weight16{0},projectable{0},nongraph{0},kernel_tests{0};};
struct Bad{bool set=false;int orbit=-1,q=-1,lin=-1,constant=-1;U support=0;std::vector<int>witnesses;};
int main(int ac,char**av){if(ac!=3){std::cerr<<"usage: verify_rm36_popcount output.json threads(12)\n";return 2;}int threads=std::atoi(av[2]);if(threads!=12)return 2;auto began=std::chrono::steady_clock::now();
  const std::array<std::uint32_t,6> reps{0,1,20,1064,1536,12640};const std::array<U,6> sizes{1,1395,54684,468720,357120,166656};
  std::array<U,20> cubic{};int ci=0;for(int a=0;a<6;++a)for(int b=a+1;b<6;++b)for(int c=b+1;c<6;++c)cubic[ci++]=monomial({a,b,c});
  std::array<U,15> quadratic{};int qi=0;for(int a=0;a<6;++a)for(int b=a+1;b<6;++b)quadratic[qi++]=monomial({a,b});
  std::array<U,6> linear{};for(int i=0;i<6;++i)linear[i]=monomial({i});
  std::vector<U> qtruth(1<<15);for(int q=1;q<(1<<15);++q){int b=__builtin_ctz(q);qtruth[q]=qtruth[q^(1<<b)]^quadratic[b];}
  std::array<U,6> ctruth{};for(int o=0;o<6;++o)for(int b=0;b<20;++b)if(reps[o]>>b&1)ctruth[o]^=cubic[b];
  auto ks=make_kernels();std::array<Stat,6> stat;Bad bad;std::mutex lock;std::atomic<int> next{0};const int jobs=6*(1<<15);std::vector<std::thread> pool;
  for(int t=0;t<threads;++t)pool.emplace_back([&]{for(;;){int start=next.fetch_add(16);if(start>=jobs)return;for(int job=start;job<std::min(jobs,start+16);++job){int o=job>>15,q=job&32767;U base=ctruth[o]^qtruth[q];for(int lin=0;lin<64;++lin){U phase=base;for(int b=0;b<6;++b)if(lin>>b&1)phase^=linear[b];int w=pop(phase);if(w!=16&&w!=48)continue;int constant=(w==48);U support=constant?~phase:phase;if(pop(support)!=16)throw std::runtime_error("weight consistency");std::uint64_t tests=0;int k=first_injective(support,ks,&tests);stat[o].weight16++;stat[o].kernel_tests+=tests;if(k>=0)stat[o].projectable++;else{stat[o].nongraph++;std::lock_guard<std::mutex>g(lock);if(!bad.set){bad.set=true;bad.orbit=o;bad.q=q;bad.lin=lin;bad.constant=constant;bad.support=support;for(const auto& kk:ks){int h=-1;if(!collision(support,kk,&h))throw std::runtime_error("certificate error");bad.witnesses.push_back(h);}}}}}}});
  for(auto&t:pool)t.join();U weighted=0,proj=0,badcount=0;for(int o=0;o<6;++o){weighted+=sizes[o]*stat[o].weight16.load();proj+=sizes[o]*stat[o].projectable.load();badcount+=sizes[o]*stat[o].nongraph.load();}const U target=232081500ULL;bool match=weighted==target;
  std::ofstream f(av[1]);if(!f)return 3;f<<"{\n\"schema\":\"a28-r019-rm36-popcount-coset-v1\",\n\"method\":\"direct ANF truth reconstruction; popcount for every linear shift; direct point/coset collisions for each kernel; no FWHT and no difference-set mask\",\n\"threads\":12,\n\"kernel_count\":651,\n\"quadratic_masks_per_orbit\":32768,\n\"weighted_weight16_count\":"<<weighted<<",\n\"macwilliams_target\":"<<target<<",\n\"macwilliams_match\":"<<(match?"true":"false")<<",\n\"weighted_projectable_count\":"<<proj<<",\n\"weighted_nongraph_count\":"<<badcount<<",\n\"orbits\":[";for(int o=0;o<6;++o){if(o)f<<",";f<<"{\"representative\":"<<reps[o]<<",\"size\":"<<sizes[o]<<",\"weight16_count\":"<<stat[o].weight16.load()<<",\"projectable_count\":"<<stat[o].projectable.load()<<",\"nongraph_count\":"<<stat[o].nongraph.load()<<",\"kernel_collision_tests\":"<<stat[o].kernel_tests.load()<<"}";}f<<"],\n\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();if(bad.set){f<<",\n\"first_nongraph\":{\"orbit\":"<<bad.orbit<<",\"quadratic_mask\":"<<bad.q<<",\"linear_mask\":"<<bad.lin<<",\"constant\":"<<bad.constant<<",\"support_mask\":\""<<hx(bad.support)<<"\",\"all_651_direct_collision_witnesses\":[";for(size_t i=0;i<bad.witnesses.size();++i){if(i)f<<",";f<<bad.witnesses[i];}f<<"]}";}f<<"\n}\n";return match?0:10;}
