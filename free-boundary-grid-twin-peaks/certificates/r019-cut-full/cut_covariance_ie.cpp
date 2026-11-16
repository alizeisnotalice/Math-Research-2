// Exact centered vertical-cut covariance under the original K=2 law.
// The independent-set recurrence is imported read-only from the r018 source.
#define main upstream_independent_peak_dp_main
#include "../../computation/independent_peak_dp.cpp"
#undef main

static void addsigned(U256& plus, U256& minus, const U192& v, int64_t coeff) {
  if (coeff >= 0) addscaled(plus, v, uint64_t(coeff));
  else addscaled(minus, v, uint64_t(-coeff));
}

// The final covariance numerators contain products of exact K=2 counts.
struct U512 { std::array<uint64_t,8> x{}; };
static void addword(U512& a, int at, uint64_t v) {
  while (v) { assert(at<8); uint64_t old=a.x[at]; a.x[at]+=v; v=(a.x[at]<old); ++at; }
}
static U512 mul256(const U256&a,const U256&b) {
  U512 z;
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) {
    __uint128_t q=(__uint128_t)a.x[i]*b.x[j];
    addword(z,i+j,uint64_t(q)); addword(z,i+j+1,uint64_t(q>>64));
  }
  return z;
}
static int cmp512(const U512&a,const U512&b) { for(int i=7;i>=0;i--) if(a.x[i]!=b.x[i]) return a.x[i]<b.x[i]?-1:1; return 0; }
static U512 sub512(U512 a,const U512&b) { uint64_t borrow=0; for(int i=0;i<8;i++){uint64_t bi=b.x[i]+borrow;bool over=bi<b.x[i];uint64_t old=a.x[i];a.x[i]-=bi;borrow=over||old<bi;}assert(!borrow);return a; }
static U512 twice(U512 a) { U512 b; for(int i=0;i<8;i++) addword(b,i,a.x[i]); for(int i=0;i<8;i++) addword(b,i,a.x[i]); return b; }
static std::string str512(U512 a) { if(a.x==std::array<uint64_t,8>{}) return "0"; std::string s; while(a.x!=std::array<uint64_t,8>{}){__uint128_t r=0;for(int i=7;i>=0;i--){__uint128_t q=(r<<64)|a.x[i];a.x[i]=uint64_t(q/10);r=q%10;}s.push_back(char('0'+int(r)));}std::reverse(s.begin(),s.end());return s; }
struct SM256 { int sign; U256 mag; };
static SM256 signmag(const U256&p,const U256&m) { int c=cmp(p,m); return {c,c<0?sub(m,p):sub(p,m)}; }
struct SM512 { int sign; U512 mag; };
static SM512 product(SM256 a,SM256 b) { return {a.sign*b.sign,mul256(a.mag,b.mag)}; }
static SM512 subtract_signed(SM512 a,SM512 b) {
  b.sign=-b.sign;
  if(!a.sign)return b; if(!b.sign)return a;
  if(a.sign==b.sign){U512 z=a.mag; for(int i=0;i<8;i++)addword(z,i,b.mag.x[i]);return {a.sign,z};}
  int c=cmp512(a.mag,b.mag); if(!c)return {0,{}};
  return c>0?SM512{a.sign,sub512(a.mag,b.mag)}:SM512{b.sign,sub512(b.mag,a.mag)};
}
static std::string strsigned512(SM512 z) { return z.sign<0?"-"+str512(z.mag):str512(z.mag); }

int main(int argc,char**argv) {
  if(argc!=2)return 2;
  n=std::atoi(argv[1]); if(n<2||n>6)return 2;
  N=n*n; alpha=(N+1)/2;
  valid.clear(); for(uint64_t r=0;r<(1ULL<<n);r++)if(!(r&(r<<1)))valid.push_back(r);
  bags.assign(alpha+1,{}); enumerate(0,0,0);
  size_t total=0;for(auto&b:bags)total+=b.size();
  std::vector<uint64_t> sets;sets.reserve(total);std::vector<size_t> start(alpha+2);
  for(int k=0;k<=alpha;k++){start[k]=sets.size();sets.insert(sets.end(),bags[k].begin(),bags[k].end());std::vector<uint64_t>().swap(bags[k]);}start[alpha+1]=sets.size();
  std::unordered_map<uint64_t,uint32_t> index;index.max_load_factor(.7);index.reserve(total*2);for(uint32_t i=0;i<sets.size();i++)index.emplace(sets[i],i);
  std::vector<U192>C(total);C[0].x[0]=1;for(int i=2;i<=N;i++)C[0]=mulsmall192(C[0],i);
  for(int k=1;k<=alpha;k++)for(size_t z=start[k];z<start[k+1];z++){
    uint64_t s=sets[z],q=s;U192 sum;while(q){uint64_t b=q&-q;q-=b;sum=add192(sum,C[index.at(s^b)]);}C[z]=divsmall192(sum,k+__builtin_popcountll(neighbors(s)));
  }

  const int cuts=n-1;
  U256 a2p,a2m;
  std::array<U256,5> lp{},lm{};
  std::array<std::array<U256,5>,5> hp{},hm{};
  for(size_t z=0;z<sets.size();z++){
    const U192&c=C[z];if(iszero(c))continue;uint64_t s=sets[z];int k=__builtin_popcountll(s);if(k<2)continue;
    bool neg=(k&1);int64_t sign=neg?-1:1;
    addsigned(a2p,a2m,c,sign*int64_t(choose(k,2)));
    std::array<int64_t,5> sum{};std::array<std::array<int64_t,5>,5> diagonal{};
    for(uint64_t q=s;q;){uint64_t b=q&-q;q-=b;int v=__builtin_ctzll(b),x=v%n;std::array<int,5> sv{};
      for(int j=0;j<cuts;j++){sv[j]=(x<=j)?-1:1;sum[j]+=sv[j];}
      for(int j=0;j<cuts;j++)for(int l=0;l<cuts;l++)diagonal[j][l]+=int64_t(sv[j])*sv[l];
    }
    for(int j=0;j<cuts;j++)addsigned(lp[j],lm[j],c,sign*int64_t(k-1)*sum[j]);
    for(int j=0;j<cuts;j++)for(int l=0;l<cuts;l++){
      int64_t h=sum[j]*sum[l]-diagonal[j][l];
      addsigned(hp[j][l],hm[j][l],c,sign*h);
    }
  }
  assert(cmp(a2p,a2m)>0);U256 A2=sub(a2p,a2m);
  for(int j=0;j<cuts;j++)for(int l=0;l<cuts;l++){assert(cmp(hp[j][l],hp[l][j])==0&&cmp(hm[j][l],hm[l][j])==0);}
  U512 denominator=twice(twice(mul256(A2,A2))); // 4 A2^2
  bool nonpositive=true;
  std::cout<<"{\"n\":"<<n<<",\"N\":"<<N<<",\"independent_sets\":"<<total<<",\"A2\":\""<<str(A2)<<"\",\"cut_marginal_numerators\":[";
  for(int j=0;j<cuts;j++){if(j)std::cout<<',';std::cout<<'"'<<signedstr(lp[j],lm[j])<<'"';}
  std::cout<<"],\"ordered_pair_cross_numerators\":[";
  for(int j=0;j<cuts;j++){if(j)std::cout<<',';std::cout<<'[';for(int l=0;l<cuts;l++){if(l)std::cout<<',';std::cout<<'"'<<signedstr(hp[j][l],hm[j][l])<<'"';}std::cout<<']';}
  std::cout<<"],\"covariance_common_denominator\":\""<<str512(denominator)<<"\",\"covariance_numerators\":[";
  for(int j=0;j<cuts;j++){if(j)std::cout<<',';std::cout<<'[';for(int l=0;l<cuts;l++){
    if(l)std::cout<<',';SM512 left=product({1,A2},signmag(hp[j][l],hm[j][l]));left.mag=twice(left.mag);
    SM512 right=product(signmag(lp[j],lm[j]),signmag(lp[l],lm[l]));SM512 cov=subtract_signed(left,right);if(cov.sign>0)nonpositive=false;
    std::cout<<'"'<<strsigned512(cov)<<'"';
  }std::cout<<']';}
  std::cout<<"],\"all_covariance_entries_nonpositive\":"<<(nonpositive?"true":"false")<<"}\n";
}
