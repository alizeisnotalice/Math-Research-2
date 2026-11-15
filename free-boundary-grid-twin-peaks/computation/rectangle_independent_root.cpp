// Exact root-band singleton count on an r x c rectangle via independent sets.
// This is deliberately a rectangle-only adaptation of the certified n<=6
// square code; it stores only independent sets, never all 2^(r*c) subsets.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

struct U192 { std::array<uint64_t,3> x{}; };
struct U256 { std::array<uint64_t,4> x{}; };
static U192 add192(U192 a,const U192&b){for(int i=0;i<3;i++){uint64_t o=a.x[i];a.x[i]+=b.x[i];if(a.x[i]<o)for(int j=i+1;j<3;j++){++a.x[j];if(a.x[j])break;}}return a;}
static U192 mulsmall192(U192 a,uint64_t q){U192 r;__uint128_t c=0;for(int i=0;i<3;i++){__uint128_t z=(__uint128_t)a.x[i]*q+c;r.x[i]=(uint64_t)z;c=z>>64;}assert(!c);return r;}
static U192 divsmall192(U192 a,uint64_t q){U192 r;__uint128_t rem=0;for(int i=2;i>=0;i--){__uint128_t z=(rem<<64)|a.x[i];r.x[i]=(uint64_t)(z/q);rem=z%q;}assert(!rem);return r;}
static bool zero(const U192&a){return !(a.x[0]|a.x[1]|a.x[2]);}
static void add256(U256&a,const U192&b){uint64_t o=a.x[0];a.x[0]+=b.x[0];uint64_t c=a.x[0]<o;for(int i=1;i<3;i++){__uint128_t z=(__uint128_t)b.x[i]+c;o=a.x[i];a.x[i]+=(uint64_t)z;c=(a.x[i]<o)+(z>>64);}o=a.x[3];a.x[3]+=c;assert(a.x[3]>=o);}
static int cmp(const U256&a,const U256&b){for(int i=3;i>=0;i--)if(a.x[i]!=b.x[i])return a.x[i]<b.x[i]?-1:1;return 0;}
static U256 sub(U256 a,const U256&b){uint64_t borrow=0;for(int i=0;i<4;i++){uint64_t bi=b.x[i]+borrow;bool o=bi<b.x[i];uint64_t old=a.x[i];a.x[i]-=bi;borrow=o||old<bi;}assert(!borrow);return a;}
static std::string dec(U256 a){if(!(a.x[0]|a.x[1]|a.x[2]|a.x[3]))return"0";std::string s;while(a.x[0]|a.x[1]|a.x[2]|a.x[3]){__uint128_t r=0;for(int i=3;i>=0;i--){__uint128_t z=(r<<64)|a.x[i];a.x[i]=(uint64_t)(z/10);r=z%10;}s.push_back(char('0'+(int)r));}std::reverse(s.begin(),s.end());return s;}
static std::string diff(const U256&p,const U256&m){return cmp(p,m)<0?"-"+dec(sub(m,p)):dec(sub(p,m));}

static int rows,cols,N,alpha;
static std::vector<uint64_t> valid;
static std::vector<std::vector<uint64_t>> bags;
static void enumerate(int row,uint64_t previous,uint64_t set){
  if(row==rows){bags[__builtin_popcountll(set)].push_back(set);return;}
  for(uint64_t r:valid)if(!(r&previous))enumerate(row+1,r,set|(r<<(row*cols)));
}
static uint64_t neighbors(uint64_t s){
  uint64_t col0=0,colLast=0;
  for(int y=0;y<rows;y++){col0|=1ULL<<(y*cols);colLast|=1ULL<<(y*cols+cols-1);}
  return ((((s&~col0)>>1)|((s&~colLast)<<1)|(s<<cols)|(s>>cols))&((1ULL<<N)-1));
}
int main(int argc,char**argv){
  if(argc!=3)return 2;rows=std::atoi(argv[1]);cols=std::atoi(argv[2]);N=rows*cols;
  if(rows<1||cols<1||rows>7||cols>7||N>35)return 2; alpha=(N+1)/2;
  for(uint64_t r=0;r<(1ULL<<cols);r++)if(!(r&(r<<1)))valid.push_back(r);
  bags.resize(alpha+1);enumerate(0,0,0);
  size_t total=0;for(auto&b:bags)total+=b.size();
  std::vector<uint64_t> sets;sets.reserve(total);std::vector<size_t> start(alpha+2);
  for(int k=0;k<=alpha;k++){start[k]=sets.size();sets.insert(sets.end(),bags[k].begin(),bags[k].end());std::vector<uint64_t>().swap(bags[k]);}start[alpha+1]=sets.size();
  std::unordered_map<uint64_t,uint32_t> index;index.max_load_factor(.70);index.reserve(total*2);
  for(uint32_t i=0;i<sets.size();i++)index.emplace(sets[i],i);
  std::vector<U192> C(total);C[0].x[0]=1;for(int i=2;i<=N;i++)C[0]=mulsmall192(C[0],i);
  for(int k=1;k<=alpha;k++)for(size_t z=start[k];z<start[k+1];z++){
    uint64_t s=sets[z],q=s;U192 sum;while(q){uint64_t bit=q&-q;q-=bit;sum=add192(sum,C[index.at(s^bit)]);}
    C[z]=divsmall192(sum,k+__builtin_popcountll(neighbors(s)));
  }
  uint64_t band=0;for(int y=0;y<rows;y++)for(int x=0;x<=(cols-1)/2;x++)band|=1ULL<<(y*cols+x);
  U256 plus,minus,a1p,a1m;
  for(size_t z=1;z<sets.size();z++){
    if(zero(C[z]))continue;int k=__builtin_popcountll(sets[z]);uint64_t w=__builtin_popcountll(sets[z]&band);
    if(k&1){for(uint64_t i=0;i<w;i++)add256(plus,C[z]);}else{for(uint64_t i=0;i<w;i++)add256(minus,C[z]);}
    if(k&1){for(int i=0;i<k;i++)add256(a1p,C[z]);}else{for(int i=0;i<k;i++)add256(a1m,C[z]);}
  }
  std::cout<<"{\"rows\":"<<rows<<",\"cols\":"<<cols<<",\"independent_sets\":"<<total
           <<",\"L\":\""<<diff(plus,minus)<<"\",\"A1\":\""<<diff(a1p,a1m)<<"\"}\n";
}
