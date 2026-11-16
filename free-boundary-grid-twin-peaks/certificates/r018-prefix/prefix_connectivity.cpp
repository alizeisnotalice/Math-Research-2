// Exact communication-class census for ordered valid prefix words.
// This is finite evidence for the swap/replace MH proposal only.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using Word = std::vector<int>;

struct DSU {
  std::vector<int> p, sz;
  explicit DSU(int n): p(n), sz(n, 1) { std::iota(p.begin(), p.end(), 0); }
  int find(int x) { while (x != p[x]) { p[x] = p[p[x]]; x = p[x]; } return x; }
  void join(int a, int b) { a=find(a); b=find(b); if(a==b) return; if(sz[a]<sz[b]) std::swap(a,b); p[b]=a; sz[a]+=sz[b]; }
};

struct Grid {
  int n, N;
  std::vector<uint32_t> adj;
  explicit Grid(int nn): n(nn), N(nn*nn), adj(N) {
    for (int y=0;y<n;y++) for (int x=0;x<n;x++) {
      int v=y*n+x;
      for (auto [dx,dy] : {std::pair{1,0}, {-1,0}, {0,1}, {0,-1}}) {
        int X=x+dx,Y=y+dy;
        if (0<=X && X<n && 0<=Y && Y<n) adj[v] |= uint32_t(1)<<(Y*n+X);
      }
    }
  }
};

static uint64_t encode(const Word& w) {
  uint64_t z=0;
  for (size_t i=0;i<w.size();i++) z |= uint64_t(w[i]) << (4*i);
  return z;
}
static Word decode(uint64_t z, int r) {
  Word w(r);
  for(int i=0;i<r;i++) w[i]=(z>>(4*i))&15;
  return w;
}
static bool valid(const Word& w, int tau, const Grid& g) {
  uint32_t seen=0;
  for(int t=1;t<=int(w.size());t++) {
    int v=w[t-1];
    if (seen & (uint32_t(1)<<v)) return false;
    bool birth=(g.adj[v]&seen)==0;
    bool required=(t==1 || t==tau);
    if (birth!=required) return false;
    seen |= uint32_t(1)<<v;
  }
  return true;
}
static std::string word_string(const Word& w, const Grid& g) {
  std::string s="[";
  for(size_t i=0;i<w.size();i++) {
    if(i) s+=",";
    s+="["+std::to_string(w[i]%g.n)+","+std::to_string(w[i]/g.n)+"]";
  }
  return s+"]";
}
static std::string feature(const Word& w, int tau, const Grid& g, const char* kind) {
  auto point=[&](int v){ return std::to_string(v%g.n)+":"+std::to_string(v/g.n); };
  if(std::string(kind)=="first") return point(w[0]);
  if(std::string(kind)=="roots") {
    if(tau>int(w.size())) return "only-first";
    std::string a=point(w[0]),b=point(w[tau-1]);
    return a<b?a+"|"+b:b+"|"+a;
  }
  if(std::string(kind)=="root-distance") {
    if(tau>int(w.size())) return "none";
    int a=w[0],b=w[tau-1];
    return std::to_string(std::abs(a%g.n-b%g.n)+std::abs(a/g.n-b/g.n));
  }
  return "";
}

struct Result {
  int n,tau,r,states,components;
  std::map<int,int> spectrum;
  std::vector<std::string> reps, invariant_candidates;
};

static Result one_case(const Grid& g, int r, int tau) {
  std::vector<uint64_t> states;
  Word w;
  auto rec = [&](auto&& self, uint32_t seen)->void {
    int t=int(w.size())+1;
    if(t>r) { states.push_back(encode(w)); return; }
    bool required=(t==1 || t==tau);
    for(int v=0;v<g.N;v++) if(!(seen&(uint32_t(1)<<v))) {
      bool birth=(g.adj[v]&seen)==0;
      if(birth!=required) continue;
      w.push_back(v); self(self, seen|(uint32_t(1)<<v)); w.pop_back();
    }
  };
  rec(rec,0);
  std::sort(states.begin(),states.end());
  std::unordered_map<uint64_t,int> where;
  where.reserve(states.size()*2+1);
  for(int i=0;i<int(states.size());i++) where.emplace(states[i],i);
  DSU dsu(int(states.size()));
  for(int a=0;a<int(states.size());a++) {
    Word x=decode(states[a],r);
    for(int i=0;i<r;i++) for(int v=0;v<g.N;v++) {
      Word y=x;
      auto it=std::find(y.begin(),y.end(),v);
      if(it==y.end()) y[i]=v; else std::swap(y[i],*it);
      auto got=where.find(encode(y));
      if(got!=where.end()) dsu.join(a,got->second);
    }
  }
  std::map<int,std::vector<int>> byroot;
  for(int i=0;i<int(states.size());i++) byroot[dsu.find(i)].push_back(i);
  Result z{g.n,tau,r,int(states.size()),int(byroot.size())};
  for(auto& [root, members]:byroot) z.spectrum[int(members.size())]++;
  if(z.components>1) {
    for(auto& [root,members]:byroot) z.reps.push_back(word_string(decode(states[members[0]],r),g));
    for(const char* key: {"first","roots","root-distance"}) {
      bool constant=true;
      std::vector<std::string> component_values;
      for(auto& [root,members]:byroot) {
        std::string f=feature(decode(states[members[0]],r),tau,g,key);
        for(int id:members) if(feature(decode(states[id],r),tau,g,key)!=f) {constant=false;break;}
        component_values.push_back(f);
        if(!constant)break;
      }
      std::sort(component_values.begin(), component_values.end());
      component_values.erase(std::unique(component_values.begin(), component_values.end()), component_values.end());
      // A global constant has no explanatory value for disconnection.
      if(constant && component_values.size()>1) z.invariant_candidates.push_back(key);
    }
  }
  return z;
}

static void json_string(std::ostream& o,const std::string& s){o<<'"';for(char c:s){if(c=='"'||c=='\\')o<<'\\';o<<c;}o<<'"';}
static void print_result(std::ostream&o,const Result&z) {
  o<<"{\"n\":"<<z.n<<",\"tau\":"<<z.tau<<",\"r\":"<<z.r<<",\"states\":"<<z.states<<",\"components\":"<<z.components<<",\"size_spectrum\":{";
  bool first=true;for(auto [a,b]:z.spectrum){if(!first)o<<",";first=false;json_string(o,std::to_string(a));o<<":"<<b;}o<<"},\"representatives\":[";
  for(size_t i=0;i<z.reps.size();i++){if(i)o<<",";json_string(o,z.reps[i]);}o<<"],\"component_constant_feature_candidates\":[";
  for(size_t i=0;i<z.invariant_candidates.size();i++){if(i)o<<",";json_string(o,z.invariant_candidates[i]);}o<<"]}";
}

int main() {
  std::vector<Result> rows;
  for(int n: {2,3}) { Grid g(n); for(int r=1;r<=g.N;r++) for(int tau=2;tau<=g.N;tau++) rows.push_back(one_case(g,r,tau)); }
  std::cout<<"{\"method\":\"exact valid ordered-prefix enumeration plus union-find over accepted swap/replace proposals\",\"proposal\":\"all (i,v) in [r]xV; invalid candidates omitted from off-diagonal graph\",\"rows\":[";
  for(size_t i=0;i<rows.size();i++){if(i)std::cout<<",";print_result(std::cout,rows[i]);}
  std::cout<<"],\"first_nontrivial_disconnected\":";
  bool found=false;
  for(const auto& z:rows) if(z.states>1 && z.components>1) {print_result(std::cout,z);found=true;break;}
  if(!found) std::cout<<"null";
  std::cout<<"}\n";
}
