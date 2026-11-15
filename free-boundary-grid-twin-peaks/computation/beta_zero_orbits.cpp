// Exact beta=0 count on n-by-n grids via D4 directed-edge orbits and
// independent-set C/X/E recurrences. Uses arbitrary precision integers.
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include <gmpxx.h>
using cpp_int=mpz_class;
using Mask=uint64_t;

static int pop(Mask x){return __builtin_popcountll(x);}
static int bit(Mask x){return __builtin_ctzll(x);}
static std::string text(const cpp_int& x){return x.get_str();}

struct Grid {
  int n,N;
  std::vector<Mask> adj,closed;
  explicit Grid(int side):n(side),N(side*side),adj(N),closed(N){
    for(int v=0;v<N;v++)for(int u=0;u<N;u++)
      if(std::abs(v%n-u%n)+std::abs(v/n-u/n)==1)adj[v]|=Mask(1)<<u;
    for(int v=0;v<N;v++)closed[v]=adj[v]|(Mask(1)<<v);
  }
  int map_d4(int v,int op) const {
    int x=v%n,y=v/n; bool reflect=op>=4; int rot=op%4;
    if(reflect)x=n-1-x;
    for(int k=0;k<rot;k++){int nx=n-1-y,ny=x;x=nx;y=ny;}
    return y*n+x;
  }
};

static cpp_int factorial(int n){cpp_int z=1;for(int i=2;i<=n;i++)z*=i;return z;}

struct Answer { cpp_int e2=0, vhigh=0; size_t sets=0; };

static Answer edge_answer(const Grid& g,int v,int x){
  Mask T=g.closed[v]|g.closed[x];
  Mask full=(Mask(1)<<g.N)-1, U=full^T;
  std::vector<Mask> masks;
  std::function<void(Mask,Mask)> rec=[&](Mask avail,Mask chosen){
    if(!avail){masks.push_back(chosen);return;}
    Mask b=avail&-avail;int z=bit(b);
    rec(avail^b,chosen);
    rec((avail^b)&~g.adj[z],chosen|b);
  };
  rec(U,0);
  std::sort(masks.begin(),masks.end(),[](Mask a,Mask b){int pa=pop(a),pb=pop(b);return pa<pb||(pa==pb&&a<b);});
  std::unordered_map<Mask,size_t> where;
  where.reserve(masks.size()*2);
  for(size_t i=0;i<masks.size();i++)where.emplace(masks[i],i);
  std::vector<cpp_int> C(masks.size()),X(masks.size()),E(masks.size());
  cpp_int fact=factorial(g.N);
  for(size_t z=0;z<masks.size();z++){
    Mask I=masks[z],cl=0;
    for(Mask q=I;q;q&=q-1){int i=bit(q);cl|=g.closed[i];}
    int csize=pop(cl),a=pop(cl|T);
    cpp_int sc=0,sx=0,se=0;
    if(!I) C[z]=fact;
    else {
      for(Mask q=I;q;q&=q-1){Mask J=I^(q&-q);sc+=C[where.at(J)];sx+=X[where.at(J)];se+=E[where.at(J)];}
      if(sc%csize!=0)throw std::runtime_error("nonintegral C recurrence");
      C[z]=sc/csize;
    }
    if(!I){sx=0;se=0;}
    if((C[z]+sx)%(a-1)!=0)throw std::runtime_error("nonintegral X recurrence");
    X[z]=(C[z]+sx)/(a-1);
    if((X[z]+se)%a!=0)throw std::runtime_error("nonintegral E recurrence");
    E[z]=(X[z]+se)/a;
    if(X[z]%g.N!=0)throw std::runtime_error("X is not divisible by N");
  }
  Answer ans;ans.sets=masks.size();
  for(size_t z=1;z<masks.size();z++){
    int k=pop(masks[z]);
    if(k&1){ans.e2+=cpp_int(k)*E[z];ans.vhigh+=cpp_int(k)*(X[z]/g.N);}
    else {ans.e2-=cpp_int(k)*E[z];ans.vhigh-=cpp_int(k)*(X[z]/g.N);}
  }
  return ans;
}

int main(int argc,char**argv){
 try {
  if(argc!=2)return 2;
  int n=std::stoi(argv[1]);
  if(n<2||n>6)return 2;
  Grid g(n);
  std::vector<std::pair<int,int>> edges;
  for(int v=0;v<g.N;v++)for(Mask q=g.adj[v];q;q&=q-1)edges.push_back({v,bit(q)});
  std::unordered_map<std::string,std::vector<std::pair<int,int>>> buckets;
  for(auto e:edges){std::string key;
    for(int op=0;op<8;op++){int a=g.map_d4(e.first,op),b=g.map_d4(e.second,op);std::string s=std::to_string(a)+":"+std::to_string(b);if(op==0||s<key)key=s;}
    buckets[key].push_back(e);
  }
  size_t orbit_weight=0;
  for(const auto& [key,orbit]:buckets)orbit_weight+=orbit.size();
  if(orbit_weight!=size_t(4*g.n*(g.n-1)))throw std::runtime_error("D4 orbit coverage mismatch");
  struct OrbitResult { std::string key; std::pair<int,int> edge; size_t multiplicity; Answer answer; cpp_int low; };
  std::vector<OrbitResult> results;
  std::unordered_map<std::string,cpp_int> low_by_key;
  cpp_int total=0;
  for(auto& [key,orbit]:buckets){
    auto e=orbit.front();Answer a=edge_answer(g,e.first,e.second);cpp_int low=a.e2-a.vhigh;
    if(low<0)throw std::runtime_error("negative beta0 edge count");
    total+=cpp_int(orbit.size())*low;
    results.push_back({key,e,orbit.size(),std::move(a),std::move(low)});
    low_by_key.emplace(key,results.back().low);
  }
  for(const auto& item:results){
    std::string reverse_key;
    for(int op=0;op<8;op++){int a=g.map_d4(item.edge.second,op),b=g.map_d4(item.edge.first,op);std::string s=std::to_string(a)+":"+std::to_string(b);if(op==0||s<reverse_key)reverse_key=s;}
    auto it=low_by_key.find(reverse_key);
    if(it==low_by_key.end()||it->second!=item.low)throw std::runtime_error("reverse-edge beta0 mismatch");
  }
  std::cout<<"{\"n\":"<<g.n<<",\"orbits\":[";bool first=true;
  for(const auto& item:results){
    if(!first)std::cout<<",";first=false;
    const auto& e=item.edge;const auto& a=item.answer;const auto& low=item.low;
    std::cout<<"{\"representative\":["<<e.first<<","<<e.second<<"],\"multiplicity\":"<<item.multiplicity
      <<",\"independent_sets\":"<<a.sets<<",\"E_K2\":\""<<text(a.e2)<<"\",\"v_high\":\""<<text(a.vhigh)<<"\",\"beta0\":\""<<text(low)<<"\"}";
  }
  std::cout<<"],\"directed_edge_coverage\":"<<orbit_weight<<",\"reverse_beta0_verified\":true,\"beta0_total\":\""<<text(total)<<"\"}\n";
  return 0;
 } catch(const std::exception& e) { std::cerr<<"beta_zero_orbits: "<<e.what()<<"\n"; return 1; }
}
