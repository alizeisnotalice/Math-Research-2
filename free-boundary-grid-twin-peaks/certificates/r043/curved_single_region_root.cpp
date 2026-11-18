// Exact one-peak root-band counts for variable-width left regions of P_n square.
// Vertices are y*n+x; a profile h gives H={ (x,y): x<h[y] }.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using U128 = unsigned __int128;
struct U256 { std::array<uint64_t,4> x{}; };

static U128 uadd(U128 a,U128 b){ assert(std::numeric_limits<U128>::max()-a>=b); return a+b; }
static U128 umul_small(U128 a,uint64_t b){ assert(b==0 || a<=std::numeric_limits<U128>::max()/b); return a*b; }
static U128 udiv_exact(U128 a,uint64_t b){ assert(b && a%b==0); return a/b; }
static std::string dec128(U128 a){ if(!a)return "0"; std::string s; while(a){s.push_back(char('0'+a%10));a/=10;} std::reverse(s.begin(),s.end());return s; }
static void add_word(U256& a,int pos,uint64_t v){ while(v){assert(pos<4);uint64_t old=a.x[pos];a.x[pos]+=v;v=a.x[pos]<old;++pos;} }
static U256 mul128(U128 a,U128 b){
  uint64_t av[2]={uint64_t(a),uint64_t(a>>64)},bv[2]={uint64_t(b),uint64_t(b>>64)};U256 r;
  for(int i=0;i<2;i++)for(int j=0;j<2;j++){U128 p=U128(av[i])*bv[j];add_word(r,i+j,uint64_t(p));add_word(r,i+j+1,uint64_t(p>>64));}
  return r;
}
static U256 mulsmall256(U256 a,uint64_t b){U128 carry=0;for(int i=0;i<4;i++){U128 p=U128(a.x[i])*b+carry;a.x[i]=uint64_t(p);carry=p>>64;}assert(!carry);return a;}
static std::string dec256(U256 a){if(!(a.x[0]|a.x[1]|a.x[2]|a.x[3]))return "0";std::string s;while(a.x[0]|a.x[1]|a.x[2]|a.x[3]){U128 rem=0;for(int i=3;i>=0;i--){U128 q=(rem<<64)|a.x[i];a.x[i]=uint64_t(q/10);rem=q%10;}s.push_back(char('0'+uint64_t(rem)));}std::reverse(s.begin(),s.end());return s;}

struct Answer { int vertices=0,allowed_vertices=0,alpha=0; uint64_t independent_sets=0; U128 rooted=0; };
static int n,N;
static std::vector<uint64_t> adj;

static uint64_t choose_small(int total,int take){
  if(take<0||take>total)return 0;take=std::min(take,total-take);U128 a=1;
  for(int i=1;i<=take;i++)a=udiv_exact(umul_small(a,uint64_t(total-take+i)),i);
  assert(a<=std::numeric_limits<uint64_t>::max());return uint64_t(a);
}
static std::vector<uint64_t> valid_row_masks(int width){
  std::vector<uint64_t> out; const uint64_t cap=uint64_t(1)<<width;
  for(uint64_t m=0;m<cap;m++)if(!(m&(m<<1)))out.push_back(m);return out;
}
static void enumerate_sets(int row,uint64_t prev,uint64_t full,const std::vector<std::vector<uint64_t>>& choices,std::vector<std::vector<uint64_t>>& bags){
  if(row==n){bags[__builtin_popcountll(full)].push_back(full);return;}
  for(uint64_t m:choices[row])if(!(m&prev))enumerate_sets(row+1,m,full|(m<<(row*n)),choices,bags);
}

static Answer solve_left(const std::string& h){
  assert(int(h.size())==n); uint64_t region=0;int vertices=0;
  std::vector<std::vector<uint64_t>> choices(n);
  for(int y=0;y<n;y++){
    int width=h[y]-'0';assert(1<=width&&width<n);vertices+=width;
    for(int x=0;x<width;x++)region|=uint64_t(1)<<(y*n+x);
    choices[y]=valid_row_masks(width);
  }
  // A connected irregular bipartite region can have an independent side larger
  // than ceil(|H|/2), so allocate all possible cardinality buckets.
  std::vector<std::vector<uint64_t>> bags(vertices+1);
  enumerate_sets(0,0,0,choices,bags);
  int alpha=0;uint64_t independent_sets=0;
  for(int k=0;k<int(bags.size());k++){if(!bags[k].empty())alpha=std::max(alpha,k);independent_sets+=bags[k].size();}
  std::vector<uint64_t> sets;sets.reserve(independent_sets);std::vector<size_t> start(alpha+2);
  for(int k=0;k<=alpha;k++){start[k]=sets.size();sets.insert(sets.end(),bags[k].begin(),bags[k].end());std::vector<uint64_t>().swap(bags[k]);}start[alpha+1]=sets.size();
  std::unordered_map<uint64_t,uint32_t> index;index.max_load_factor(.70);index.reserve(sets.size()*2);
  for(uint32_t i=0;i<sets.size();i++)index.emplace(sets[i],i);
  std::vector<U128> C(sets.size());C[0]=1;for(int i=2;i<=vertices;i++)C[0]=umul_small(C[0],i);
  for(int k=1;k<=alpha;k++)for(size_t z=start[k];z<start[k+1];z++){
    uint64_t s=sets[z],q=s,neigh=0;while(q){uint64_t b=q&-q;q-=b;neigh|=adj[__builtin_ctzll(b)]&region;}
    const int denom=__builtin_popcountll(s|neigh);U128 sum=0;q=s;
    while(q){uint64_t b=q&-q;q-=b;auto it=index.find(s^b);assert(it!=index.end());sum=uadd(sum,C[it->second]);}
    C[z]=udiv_exact(sum,denom);
  }
  const int band=(n+3)/4;uint64_t allowed=0;
  for(int y=0;y<n;y++)for(int x=0;x<band;x++){int v=y*n+x;uint64_t b=uint64_t(1)<<v;if((region&b)&&!(adj[v]&~region))allowed|=b;}
  U128 plus=0,minus=0;
  for(size_t z=1;z<sets.size();z++){
    const uint64_t s=sets[z];const uint64_t w=__builtin_popcountll(s&allowed);if(!w)continue;
    U128 term=umul_small(C[z],w);if(__builtin_popcountll(s)&1)plus=uadd(plus,term);else minus=uadd(minus,term);
  }
  assert(plus>=minus);
  return {vertices,int(__builtin_popcountll(allowed)),alpha,independent_sets,plus-minus};
}
static std::string complement_profile(const std::string& h){std::string r=h;for(char&c:r)c=char('0'+n-(c-'0'));return r;}

int main(int argc,char**argv){
  if(argc!=4){std::cerr<<"usage: N OUTPUT PROFILE\\n";return 2;}
  n=std::atoi(argv[1]);N=n*n;if((n!=5&&n!=7)||N>=64)return 2;const std::string h=argv[3];
  const int low=n==5?2:3,high=n==5?3:4;
  if(int(h.size())!=n||std::filesystem::exists(argv[2]))return 2;
  for(char c:h)if(c<'0'+low||c>'0'+high)return 2;
  adj.assign(N,0);for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(std::abs(v%n-u%n)+std::abs(v/n-u/n)==1)adj[v]|=uint64_t(1)<<u;
  Answer left=solve_left(h);const std::string hc=complement_profile(h);Answer right=solve_left(hc);
  uint64_t shuffle=choose_small(N,left.vertices);U256 event=mulsmall256(mul128(left.rooted,right.rooted),shuffle);
  std::ofstream out(argv[2]);if(!out)return 5;
  out<<"{\"schema\":\"curved-single-region-root-ie-v1\",\"n\":"<<n<<",\"curve\":\""<<h<<"\",\"complement_curve\":\""<<hc<<"\","
     <<"\"left_vertices\":"<<left.vertices<<",\"right_vertices\":"<<right.vertices<<",\"left_allowed_vertices\":"<<left.allowed_vertices<<",\"right_allowed_vertices\":"<<right.allowed_vertices<<","
     <<"\"left_independent_sets\":"<<left.independent_sets<<",\"right_independent_sets\":"<<right.independent_sets<<",\"left_alpha\":"<<left.alpha<<",\"right_alpha\":"<<right.alpha<<","
     <<"\"left_rooted_count\":\""<<dec128(left.rooted)<<"\",\"right_rooted_count\":\""<<dec128(right.rooted)<<"\",\"shuffle\":\""<<shuffle<<"\",\"count\":\""<<dec256(event)<<"\"}\n";
  return 0;
}
