// Complete orbit-compressed exact audit of RM(3,6) weight-16 supports.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using U=std::uint64_t;
static int pc(U x){int n=0;while(x){n+=x&1;x>>=1;}return n;}
static int triple_index(int a,int b,int c,const std::vector<std::array<int,3>>& T){for(int i=0;i<(int)T.size();++i)if(T[i]==std::array<int,3>{a,b,c})return i;return -1;}
static U truth_monomial(const std::vector<int>& vars){U m=0;for(int x=0;x<64;++x){int v=1;for(int z:vars)v&=(x>>z)&1;if(v)m|=1ULL<<x;}return m;}
static void fwht(std::array<int,64>& a){for(int h=1;h<64;h<<=1)for(int p=0;p<64;p+=2*h)for(int i=p;i<p+h;++i){int u=a[i],v=a[i+h];a[i]=u+v;a[i+h]=u-v;}}
static std::string hex(U x){std::ostringstream o;o<<"0x"<<std::hex<<x;return o.str();}

struct Kernel{int a,b,c;U mask;};
static std::vector<Kernel> kernels(){std::set<int> seen;std::vector<Kernel> out;for(int a=1;a<64;++a)for(int b=a+1;b<64;++b){int c=a^b;std::array<int,3> q{a,b,c};std::sort(q.begin(),q.end());int key=q[0]|(q[1]<<6)|(q[2]<<12);if(seen.insert(key).second)out.push_back({q[0],q[1],q[2],(1ULL<<q[0])|(1ULL<<q[1])|(1ULL<<q[2])});}if(out.size()!=651)throw std::runtime_error("kernel count");return out;}
static U difference_mask(U support){U d=1;for(int x=0;x<64;++x)if((support>>x)&1)for(int y=0;y<64;++y)if((support>>y)&1)d|=1ULL<<(x^y);return d;}
static bool direct_collision(U support,const Kernel& k){for(int x=0;x<64;++x)if((support>>x)&1)for(int h:{k.a,k.b,k.c})if((support>>(x^h))&1)return true;return false;}
static int projecting_kernel(U diff,const std::vector<Kernel>& K,std::uint64_t* scanned=nullptr){for(int i=0;i<(int)K.size();++i){if(scanned)++*scanned;if((diff&K[i].mask)==0)return i;}return -1;}
static void check_controls(const std::vector<Kernel>& K){
  U positive=(1ULL<<16)-1, negative=0;U state=0x9e3779b97f4a7c15ULL;
  for(int attempt=0;attempt<10000 && !negative;++attempt){U s=0;while(pc(s)<16){state^=state<<7;state^=state>>9;state^=state<<8;s|=1ULL<<(state&63);}if(projecting_kernel(difference_mask(s),K)<0)negative=s;}
  if(!negative)throw std::runtime_error("no deterministic negative control");
  for(U s:{positive,negative}){U d=difference_mask(s);for(const auto& k:K)if(((d&k.mask)!=0)!=direct_collision(s,k))throw std::runtime_error("difference/direct collision mismatch");}
  if(projecting_kernel(difference_mask(positive),K)<0 || projecting_kernel(difference_mask(negative),K)>=0)throw std::runtime_error("control polarity");
}

struct Orbit{std::uint32_t rep,size;};
struct Maps{std::array<std::uint32_t,20> cyc,swp,tr;};
static std::uint32_t act(std::uint32_t m,const std::array<std::uint32_t,20>& map){std::uint32_t o=0;while(m){int b=__builtin_ctz(m);o^=map[b];m&=m-1;}return o;}
static Maps make_maps(const std::vector<std::array<int,3>>& T){Maps M{};for(int i=0;i<20;++i){auto t=T[i];auto put=[&](std::array<int,3> q){std::sort(q.begin(),q.end());int j=triple_index(q[0],q[1],q[2],T);if(j>=0)return 1U<<j;return 0U;};auto c=t;for(int&v:c)v=(v+1)%6;M.cyc[i]=put(c);auto s=t;for(int&v:s){if(v==0)v=1;else if(v==1)v=0;}M.swp[i]=put(s);M.tr[i]=1U<<i;if(std::find(t.begin(),t.end(),0)!=t.end()){auto q=t;for(int&v:q)if(v==0)v=1;M.tr[i]^=put(q);}}return M;}
static std::vector<Orbit> orbits(const Maps& M){const int ALL=1<<20;std::vector<unsigned char> seen(ALL);std::vector<Orbit> out;for(int seed=0;seed<ALL;++seed)if(!seen[seed]){std::vector<std::uint32_t> q{(std::uint32_t)seed};seen[seed]=1;for(size_t h=0;h<q.size();++h)for(const auto* g:{&M.cyc,&M.swp,&M.tr}){auto y=act(q[h],*g);if(!seen[y]){seen[y]=1;q.push_back(y);}}out.push_back({(std::uint32_t)seed,(std::uint32_t)q.size()});}U sum=0;for(auto o:out)sum+=o.size;if(sum!=ALL)throw std::runtime_error("orbit coverage");return out;}
static void validate_maps(const Maps& M){for(std::uint32_t x=0;x<(1U<<20);++x){if(act(act(x,M.swp),M.swp)!=x||act(act(x,M.tr),M.tr)!=x)throw std::runtime_error("involution");auto y=x;for(int k=0;k<6;++k)y=act(y,M.cyc);if(y!=x)throw std::runtime_error("cycle order");}}

struct Stat{std::atomic<U> candidates{0},projectable{0},nongraph{0},kernel_scans{0};};
struct Bad{bool has=false;int orbit=-1,quad=-1,lin=-1,constant=-1;U support=0,diff=0;std::vector<int> witnesses;};
static std::string anf(std::uint32_t cub,std::uint32_t quad,int lin,int c,const std::vector<std::array<int,3>>& T,const std::vector<std::array<int,2>>& Q){std::vector<std::string> z;for(int i=0;i<20;++i)if(cub>>i&1)z.push_back("x"+std::to_string(T[i][0])+"x"+std::to_string(T[i][1])+"x"+std::to_string(T[i][2]));for(int i=0;i<15;++i)if(quad>>i&1)z.push_back("x"+std::to_string(Q[i][0])+"x"+std::to_string(Q[i][1]));for(int i=0;i<6;++i)if(lin>>i&1)z.push_back("x"+std::to_string(i));if(c)z.push_back("1");std::ostringstream o;for(size_t i=0;i<z.size();++i){if(i)o<<"+";o<<z[i];}return o.str();}

template<class F> static void parallel(int jobs,int threads,F fn){std::atomic<int> next{0};std::vector<std::thread> p;for(int t=0;t<threads;++t)p.emplace_back([&]{for(;;){int k=next.fetch_add(16);if(k>=jobs)return;for(int j=k;j<std::min(jobs,k+16);++j)fn(j);}});for(auto&x:p)x.join();}

int main(int ac,char**av){if(ac!=4){std::cerr<<"usage: audit_rm36_weight16 orbit|full threads output.json\n";return 2;}std::string mode=av[1];int threads=std::atoi(av[2]);if((mode!="orbit"&&mode!="full")||threads<1||threads>15)return 2;
  std::vector<std::array<int,3>> T;std::vector<std::array<int,2>> Q;for(int a=0;a<6;++a)for(int b=a+1;b<6;++b){Q.push_back({a,b});for(int c=b+1;c<6;++c)if(a==0||true){} }for(int a=0;a<6;++a)for(int b=a+1;b<6;++b)for(int c=b+1;c<6;++c)T.push_back({a,b,c});
  if(T.size()!=20||Q.size()!=15)throw std::runtime_error("monomial counts");Maps M=make_maps(T);validate_maps(M);auto O=orbits(M);auto K=kernels();check_controls(K);U osum=0;for(auto o:O)osum+=o.size;auto began=std::chrono::steady_clock::now();
  if(mode=="orbit"){std::array<U,20> ct{};std::array<U,15> qt{};for(int i=0;i<20;++i)ct[i]=truth_monomial({T[i][0],T[i][1],T[i][2]});for(int i=0;i<15;++i)qt[i]=truth_monomial({Q[i][0],Q[i][1]});std::ofstream f(av[3]);f<<"{\n\"schema\":\"a28-r019-rm36-orbits-v2\",\n\"variable_count\":6,\n\"truth_table_size\":64,\n\"orbit_count\":"<<O.size()<<",\n\"orbit_size_sum\":"<<osum<<",\n\"cubic_basis_order\":[";for(int i=0;i<20;++i){if(i)f<<",";f<<"\"x"<<T[i][0]<<"x"<<T[i][1]<<"x"<<T[i][2]<<"\"";}f<<"],\n\"quadratic_basis_order\":[";for(int i=0;i<15;++i){if(i)f<<",";f<<"\"x"<<Q[i][0]<<"x"<<Q[i][1]<<"\"";}f<<"],\n\"quadratic_truth_hex64\":[";for(int i=0;i<15;++i){if(i)f<<",";f<<"\""<<hex(qt[i])<<"\"";}f<<"],\n\"orbit_representatives\":[";for(size_t i=0;i<O.size();++i){if(i)f<<",";U truth=0;for(int b=0;b<20;++b)if(O[i].rep>>b&1)truth^=ct[b];f<<"{\"cubic_coefficient_mask\":"<<O[i].rep<<",\"cubic_truth_hex64\":\""<<hex(truth)<<"\",\"orbit_size\":"<<O[i].size<<"}";}f<<"],\n\"generator_validation\":\"cycle^6=swap^2=transvection^2=id on all 2^20 masks\",\n\"kernel_count\":651,\n\"control_difference_vs_direct_collision\":true,\n\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<"\n}\n";return 0;}
  std::array<U,20> ct{};std::array<U,15> qt{};std::array<U,6> lt{};for(int i=0;i<20;++i)ct[i]=truth_monomial({T[i][0],T[i][1],T[i][2]});for(int i=0;i<15;++i)qt[i]=truth_monomial({Q[i][0],Q[i][1]});for(int i=0;i<6;++i)lt[i]=truth_monomial({i});
  std::vector<U> quad(1<<15);for(int s=1;s<(1<<15);++s){int b=__builtin_ctz(s);quad[s]=quad[s^(1<<b)]^qt[b];}std::vector<U> cubic(O.size());for(size_t o=0;o<O.size();++o)for(int b=0;b<20;++b)if(O[o].rep>>b&1)cubic[o]^=ct[b];
  std::unique_ptr<Stat[]> stats(new Stat[O.size()]);Bad bad;std::mutex badlock;int jobs=(int)O.size()*(1<<15);
  parallel(jobs,threads,[&](int job){int oi=job>>15,q=job&((1<<15)-1);U base=cubic[oi]^quad[q];std::array<int,64>w;for(int x=0;x<64;++x)w[x]=(base>>x&1)?-1:1;fwht(w);for(int lin=0;lin<64;++lin)if(std::abs(w[lin])==32){int constant=(w[lin]<0);U support=base^(constant?~0ULL:0ULL);for(int b=0;b<6;++b)if(lin>>b&1)support^=lt[b];if(pc(support)!=16)throw std::runtime_error("weight relation");U diff=difference_mask(support);U scans=0;int ki=projecting_kernel(diff,K,&scans);stats[oi].candidates++;stats[oi].kernel_scans+=scans;if(ki>=0)stats[oi].projectable++;else{stats[oi].nongraph++;std::lock_guard<std::mutex> lock(badlock);if(!bad.has){bad.has=true;bad.orbit=oi;bad.quad=q;bad.lin=lin;bad.constant=constant;bad.support=support;bad.diff=diff;for(const auto& k:K){int witness=-1;for(int h:{k.a,k.b,k.c})if(diff>>h&1){witness=h;break;}bad.witnesses.push_back(witness);}}}}});
  U weighted=0,weightedproj=0,weightedbad=0,totalcand=0;for(size_t i=0;i<O.size();++i){weighted+=stats[i].candidates.load()*O[i].size;weightedproj+=stats[i].projectable.load()*O[i].size;weightedbad+=stats[i].nongraph.load()*O[i].size;totalcand+=stats[i].candidates.load();}const U macwilliams_target=232081500ULL;bool macwilliams_match=(weighted==macwilliams_target);
  std::ofstream f(av[3]);f<<"{\n\"schema\":\"a28-r019-rm36-weight16-v1\",\n\"orbit_count\":"<<O.size()<<",\n\"orbit_size_sum\":"<<osum<<",\n\"quadratic_masks_per_orbit\":32768,\n\"kernel_count\":651,\n\"total_weight16_count_weighted_by_orbit\":"<<weighted<<",\n\"macwilliams_independent_target\":"<<macwilliams_target<<",\n\"macwilliams_count_match\":"<<(macwilliams_match?"true":"false")<<",\n\"total_projectable_weight16_count_weighted_by_orbit\":"<<weightedproj<<",\n\"total_nongraph_weight16_count_weighted_by_orbit\":"<<weightedbad<<",\n\"orbits\":[";for(size_t i=0;i<O.size();++i){if(i)f<<",";f<<"{\"representative\":"<<O[i].rep<<",\"size\":"<<O[i].size<<",\"quadratic_coverage\":32768,\"weight16_representative_count\":"<<stats[i].candidates.load()<<",\"projectable_representative_count\":"<<stats[i].projectable.load()<<",\"nongraph_representative_count\":"<<stats[i].nongraph.load()<<",\"kernel_tests_examined\":"<<stats[i].kernel_scans.load()<<"}";}f<<"],\n\"all_weight16_projectable\":"<<(weightedbad==0?"true":"false")<<",\n\"generator_justification\":\"permutations from cycle and swap; conjugates of x0->x0+x1 give all elementary transvections; these generate GL(6,F2) by row reduction\",\n\"top_degree_quotient_justification\":\"linear substitutions act on cubic leading forms; generated lower-degree terms are absorbed by complete degree<=2 enumeration\",\n\"difference_mask_is_bitwise_or_including_zero\":true,\n\"control_difference_vs_direct_collision\":true,\n\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
  if(bad.has){f<<",\n\"first_nongraph\":{\"orbit\":"<<bad.orbit<<",\"cubic_representative\":"<<O[bad.orbit].rep<<",\"quadratic_mask\":"<<bad.quad<<",\"linear_mask\":"<<bad.lin<<",\"constant\":"<<bad.constant<<",\"support_mask\":\""<<hex(bad.support)<<"\",\"difference_mask\":\""<<hex(bad.diff)<<"\",\"anf\":\""<<anf(O[bad.orbit].rep,bad.quad,bad.lin,bad.constant,T,Q)<<"\",\"all_651_kernel_witness_differences\":[";for(size_t i=0;i<bad.witnesses.size();++i){if(i)f<<",";f<<bad.witnesses[i];}f<<"]}";}f<<"\n}\n";return macwilliams_match?0:10;}
