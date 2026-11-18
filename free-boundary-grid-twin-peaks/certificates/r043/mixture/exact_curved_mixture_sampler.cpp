// Exact finite sampler for a fixed safe-root curved-cut event E_h.
// Word order is descending label order. PRNG samples are numerical evidence only.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using U128=unsigned __int128;
struct U256{std::array<uint64_t,4>x{};};
static U128 add(U128 a,U128 b){if(std::numeric_limits<U128>::max()-a<b)throw std::overflow_error("u128 add");return a+b;}
static U128 mulsmall(U128 a,uint64_t b){if(b&&a>std::numeric_limits<U128>::max()/b)throw std::overflow_error("u128 mul");return a*b;}
static U128 divexact(U128 a,uint64_t b){if(!b||a%b)throw std::runtime_error("nonintegral completion recurrence");return a/b;}
static std::string dec128(U128 x){if(!x)return"0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
static void addword(U256&a,int p,uint64_t v){while(v){if(p>=4)throw std::overflow_error("u256 product");uint64_t old=a.x[p];a.x[p]+=v;v=a.x[p]<old;++p;}}
static U256 mul128(U128 a,U128 b){uint64_t A[2]={uint64_t(a),uint64_t(a>>64)},B[2]={uint64_t(b),uint64_t(b>>64)};U256 r;for(int i=0;i<2;i++)for(int j=0;j<2;j++){U128 z=U128(A[i])*B[j];addword(r,i+j,uint64_t(z));addword(r,i+j+1,uint64_t(z>>64));}return r;}
static U256 mulsmall256(U256 a,uint64_t b){U128 carry=0;for(int i=0;i<4;i++){U128 z=U128(a.x[i])*b+carry;a.x[i]=uint64_t(z);carry=z>>64;}if(carry)throw std::overflow_error("u256 mul");return a;}
static std::string dec256(U256 a){if(!(a.x[0]|a.x[1]|a.x[2]|a.x[3]))return"0";std::string s;while(a.x[0]|a.x[1]|a.x[2]|a.x[3]){U128 rem=0;for(int i=3;i>=0;i--){U128 z=(rem<<64)|a.x[i];a.x[i]=uint64_t(z/10);rem=z%10;}s.push_back(char('0'+uint64_t(rem)));}std::reverse(s.begin(),s.end());return s;}
static uint64_t choose_small(int N,int k){k=std::min(k,N-k);U128 x=1;for(int i=1;i<=k;i++)x=divexact(mulsmall(x,uint64_t(N-k+i)),i);if(x>std::numeric_limits<uint64_t>::max())throw std::overflow_error("binom");return uint64_t(x);}

struct RNG{
 uint64_t s;
 explicit RNG(uint64_t seed):s(seed?seed:0x9e3779b97f4a7c15ULL){}
 uint64_t next64(){uint64_t z=(s+=0x9e3779b97f4a7c15ULL);z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);}
 U128 below(U128 b){if(!b)throw std::runtime_error("randbelow(0)");const U128 cutoff=(-b)%b;for(;;){uint64_t hi=next64();uint64_t lo=next64();U128 x=(U128(hi)<<64)|lo;if(x>=cutoff)return x%b;}}
};

static int n,N; static std::vector<uint64_t> adj; static size_t temp_set_cap,temp_set_count;
static uint64_t neighbor(uint64_t S,uint64_t within){uint64_t out=0;while(S){uint64_t b=S&-S;S-=b;out|=adj[__builtin_ctzll(b)];}return out&within;}
static bool connected(uint64_t S){if(!S)return false;uint64_t todo=S,seen=todo&-todo,front=seen;while(front){uint64_t b=front&-front;front-=b;uint64_t more=adj[__builtin_ctzll(b)]&S&~seen;seen|=more;front|=more;}return seen==S;}
static void enumerate_independent(uint64_t avail,uint64_t picked,std::vector<std::vector<uint64_t>>&bags){
 if(!avail){if(++temp_set_count>temp_set_cap)throw std::runtime_error("temporary independent-set cap exceeded");bags[__builtin_popcountll(picked)].push_back(picked);return;}
 uint64_t v=avail&-avail;enumerate_independent(avail^v,picked,bags);int i=__builtin_ctzll(v);enumerate_independent(avail&~v&~adj[i],picked|v,bags);
}

class CompletionCounter{
 uint64_t H,allowed; size_t memo_cap; std::unordered_map<uint64_t,U128> memo;
 public:
 CompletionCounter(uint64_t h,uint64_t a,size_t cap):H(h),allowed(a),memo_cap(cap){memo.max_load_factor(.70);memo.reserve(std::min<size_t>(cap,4096));}
 U128 h(uint64_t S){
  if(!S||S&~H||!connected(S))throw std::runtime_error("invalid nonempty connected prefix");
  auto found=memo.find(S);if(found!=memo.end())return found->second;
  if(memo.size()>=memo_cap)throw std::runtime_error("prefix completion cache cap exceeded");
  uint64_t R=H&~S;if(!R)return memo.emplace(S,1).first->second;
  uint64_t F=R&~neighbor(S,H); // vertices needing a higher neighbor from the remaining word
  std::vector<std::vector<uint64_t>> bags(__builtin_popcountll(F)+1);temp_set_count=0;enumerate_independent(F,0,bags);
  int alpha=0;size_t total=0;for(int k=0;k<int(bags.size());k++)if(!bags[k].empty()){alpha=k;total+=bags[k].size();}
  std::vector<uint64_t> sets;sets.reserve(total);std::vector<size_t> start(alpha+2);
  for(int k=0;k<=alpha;k++){start[k]=sets.size();sets.insert(sets.end(),bags[k].begin(),bags[k].end());std::vector<uint64_t>().swap(bags[k]);}start[alpha+1]=sets.size();
  std::unordered_map<uint64_t,uint32_t> idx;idx.max_load_factor(.70);idx.reserve(sets.size()*2);for(uint32_t i=0;i<sets.size();i++)idx.emplace(sets[i],i);
  std::vector<U128>C(sets.size());C[0]=1;for(int i=2;i<=__builtin_popcountll(R);i++)C[0]=mulsmall(C[0],i);
  for(int k=1;k<=alpha;k++)for(size_t z=start[k];z<start[k+1];z++){
   uint64_t I=sets[z],q=I;U128 sum=0;while(q){uint64_t b=q&-q;q-=b;auto it=idx.find(I^b);if(it==idx.end())throw std::runtime_error("missing IE subset");sum=add(sum,C[it->second]);}
   C[z]=divexact(sum,__builtin_popcountll(I|neighbor(I,R)));
  }
  U128 plus=0,minus=0;for(size_t z=0;z<sets.size();z++){if(__builtin_popcountll(sets[z])&1)minus=add(minus,C[z]);else plus=add(plus,C[z]);}
  if(plus<minus)throw std::runtime_error("negative completion count");U128 ans=plus-minus;memo.emplace(S,ans);return ans;
 }
 std::vector<int> roots()const{std::vector<int>r;for(uint64_t q=allowed;q;q&=q-1)r.push_back(__builtin_ctzll(q));return r;}
 uint64_t region()const{return H;} size_t cache_size()const{return memo.size();}
};

static uint64_t region_for(const std::string&h){uint64_t H=0;for(int y=0;y<n;y++)for(int x=0;x<h[y]-'0';x++)H|=uint64_t(1)<<(y*n+x);return H;}
static std::string complement(const std::string&h){std::string r=h;for(char&c:r)c=char('0'+n-(c-'0'));return r;}
static uint64_t safe_allowed(uint64_t H){uint64_t A=0;int band=(n+3)/4;for(int y=0;y<n;y++)for(int x=0;x<band;x++){int v=y*n+x;uint64_t b=uint64_t(1)<<v;if((H&b)&&!(adj[v]&~H))A|=b;}return A;}
static std::pair<int,U128> choose_root(CompletionCounter&c,RNG&rng){std::vector<std::pair<int,U128>>w;U128 total=0;for(int v:c.roots()){U128 q=c.h(uint64_t(1)<<v);w.push_back({v,q});total=add(total,q);}if(!total)throw std::runtime_error("empty root completion mass");U128 u=rng.below(total);for(auto [v,q]:w){if(u<q)return{v,total};u-=q;}throw std::runtime_error("root draw");}
static std::vector<int> sample_side(CompletionCounter&c,RNG&rng){auto [root,mass]=choose_root(c,rng);std::vector<int>word{root};uint64_t S=uint64_t(1)<<root;while(S!=c.region()){
 uint64_t B=neighbor(S,c.region())&~S;U128 total=0;std::vector<std::pair<int,U128>>w;while(B){uint64_t b=B&-B;B-=b;int v=__builtin_ctzll(b);U128 q=c.h(S|b);w.push_back({v,q});total=add(total,q);}if(!total||total!=c.h(S))throw std::runtime_error("transition partition failure");U128 u=rng.below(total);int chosen=-1;for(auto [v,q]:w){if(u<q){chosen=v;break;}u-=q;}if(chosen<0)throw std::runtime_error("transition draw");S|=uint64_t(1)<<chosen;word.push_back(chosen);
 }return word;}
static std::vector<int> shuffle_words(const std::vector<int>&L,const std::vector<int>&R,RNG&rng){int total=L.size()+R.size();std::vector<int>pos(total);for(int i=0;i<total;i++)pos[i]=i;for(int i=total-1;i>0;i--)std::swap(pos[i],pos[size_t(rng.below(i+1))]);std::sort(pos.begin(),pos.begin()+L.size());std::vector<char>isleft(total);for(int i=0;i<int(L.size());i++)isleft[pos[i]]=1;std::vector<int>w;int a=0,b=0;for(int i=0;i<total;i++)w.push_back(isleft[i]?L[a++]:R[b++]);return w;}

int main(int argc,char**argv){
 if(argc<4){std::cerr<<"usage: N PROFILE OUTPUT [--prefix-mask HEX] [--region-only] [--samples K] [--seed U64] [--cache-mib M]\\n";return 2;}
 n=std::atoi(argv[1]);N=n*n;const std::string profile=argv[2];const char*outname=argv[3];
 int lo=0,hi=0;if(n==3){lo=1;hi=2;}else if(n==5){lo=2;hi=3;}else if(n==7){lo=3;hi=4;}else return 2;if(int(profile.size())!=n||std::filesystem::exists(outname))return 4;for(char c:profile)if(c<'0'+lo||c>'0'+hi)return 2;
 uint64_t prefix=0,seed=1;int samples=0,cache_mib=128;bool query=false,region_only=false;
 for(int i=4;i<argc;i++){std::string a=argv[i];if(a=="--prefix-mask"&&++i<argc){prefix=std::stoull(argv[i],nullptr,16);query=true;}else if(a=="--region-only")region_only=true;else if(a=="--samples"&&++i<argc)samples=std::stoi(argv[i]);else if(a=="--seed"&&++i<argc)seed=std::stoull(argv[i]);else if(a=="--cache-mib"&&++i<argc)cache_mib=std::stoi(argv[i]);else return 2;}
 if(cache_mib<8||cache_mib>128||samples<0||(region_only&&samples))return 2;temp_set_cap=size_t(cache_mib)*(size_t(1)<<20)/64;size_t memo_cap=size_t(cache_mib)*(size_t(1)<<20)/64;
 adj.assign(N,0);for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(std::abs(v%n-u%n)+std::abs(v/n-u/n)==1)adj[v]|=uint64_t(1)<<u;
 uint64_t L=region_for(profile),R=region_for(complement(profile));CompletionCounter left(L,safe_allowed(L),memo_cap);
 U128 aL=0;for(int v:left.roots())aL=add(aL,left.h(uint64_t(1)<<v));U128 hp=0;if(query){if(!prefix||prefix&~L||!connected(prefix))return 2;hp=left.h(prefix);}
 U128 aR=0;U256 event{};std::vector<std::vector<int>> words;size_t right_cache=0;int right_roots=0;
 if(!region_only){CompletionCounter right(R,safe_allowed(R),memo_cap);for(int v:right.roots())aR=add(aR,right.h(uint64_t(1)<<v));right_roots=right.roots().size();if(samples&&(!aL||!aR))throw std::runtime_error("zero event mass");if(aL&&aR)event=mulsmall256(mul128(aL,aR),choose_small(N,__builtin_popcountll(L)));RNG rng(seed);for(int t=0;t<samples;t++){auto lw=sample_side(left,rng),rw=sample_side(right,rng);for(int&v:rw)v=(v/n)*n+(n-1-v%n);words.push_back(shuffle_words(lw,rw,rng));}right_cache=right.cache_size();}
 std::ofstream out(outname);if(!out)return 5;out<<"{\"schema\":\"exact-curved-mixture-sampler-v1\",\"n\":"<<n<<",\"curve\":\""<<profile<<"\",\"region_only\":"<<(region_only?"true":"false")<<",\"complement_curve\":\""<<complement(profile)<<"\",\"left_vertices\":"<<__builtin_popcountll(L)<<",\"right_vertices\":"<<__builtin_popcountll(R)<<",\"left_allowed_roots\":"<<left.roots().size()<<",\"right_allowed_roots\":"<<right_roots<<",\"left_rooted_count\":\""<<dec128(aL)<<"\",\"right_rooted_count\":\""<<dec128(aR)<<"\",\"count\":\""<<dec256(event)<<"\",\"cache_mib\":"<<cache_mib<<",\"left_cache_entries\":"<<left.cache_size()<<",\"right_cache_entries\":"<<right_cache;if(query)out<<",\"prefix_mask\":\""<<std::hex<<prefix<<std::dec<<"\",\"prefix_completion\":\""<<dec128(hp)<<"\"";out<<",\"samples\":[";for(size_t i=0;i<words.size();i++){if(i)out<<',';out<<'[';for(size_t j=0;j<words[i].size();j++){if(j)out<<',';out<<words[i][j];}out<<']';}out<<"]}\n";
}
