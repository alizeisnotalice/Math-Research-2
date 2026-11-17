// Modular, layer-rolling selected-prefix IE.  One invocation is one prime
// residue checkpoint.  n=8 requires explicit --enable-n8 and is never the
// default; this file is first validated through n=7.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using Mask=uint64_t; using u64=uint64_t; using u128=__uint128_t;
static constexpr std::array<u64,4> PRIMES={
  9223372036854775783ULL, 9223372036854775643ULL,
  9223372036854775549ULL, 9223372036854775507ULL};
static constexpr u64 PAYLOAD_LIMIT_BYTES=7ULL<<30;
static int n,N,rows,alpha; static Mask S; static std::array<Mask,64> closed{};
static std::vector<Mask> valid; static u64 memo[7][256][25];
static u64 modp; static std::array<u64,65> inv{};
static int pc(Mask x){return __builtin_popcountll(x);} static int bit(Mask x){return __builtin_ctzll(x);}
static u64 addmod(u64 a,u64 b){u64 z=a+b;return z>=modp?z-modp:z;}
static u64 mulmod(u64 a,u64 b){return u64((u128(a)*b)%modp);}
static u64 powmod(u64 a,u64 e){u64 r=1;while(e){if(e&1)r=mulmod(r,a);a=mulmod(a,a);e>>=1;}return r;}

struct Hash {
  static constexpr Mask EMPTY=~Mask(0); std::vector<Mask> key; std::vector<uint32_t> at; size_t cap=0;
  static u64 mix(u64 x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);}
  static size_t slots_for(size_t want){if(!want)return 0;size_t z=1;while(u128(5)*z<u128(8)*want)z<<=1;return z;}
  static u64 bytes_for(size_t want){return u64(slots_for(want))*12ULL;}
  explicit Hash(size_t want=0){if(!want)return;cap=slots_for(want);key.assign(cap,EMPTY);at.resize(cap);}
  u64 bytes()const{return u64(key.size())*8ULL+u64(at.size())*4ULL;}
  void insert(Mask x,uint32_t i){size_t p=mix(x)&(cap-1);while(key[p]!=EMPTY){if(key[p]==x)throw std::runtime_error("duplicate independent mask");p=(p+1)&(cap-1);}key[p]=x;at[p]=i;}
  uint32_t find(Mask x)const{size_t p=mix(x)&(cap-1);while(key[p]!=EMPTY){if(key[p]==x)return at[p];p=(p+1)&(cap-1);}throw std::runtime_error("missing previous-layer mask");}
};

static u64 suffix(int row,Mask previous,int rem){
  if(rem<0||rem>alpha||row<0||row>rows||previous>255)throw std::runtime_error("suffix bounds");
  return memo[row][previous][rem];
}
static void precompute_suffix(){
  for(int previous=0;previous<(1<<n);previous++)for(int rem=0;rem<=alpha;rem++)memo[rows][previous][rem]=(rem==0);
  for(int row=rows-1;row>=0;row--)for(int previous=0;previous<(1<<n);previous++)for(int rem=0;rem<=alpha;rem++){u64 total=0;for(Mask r:valid)if(!(r&Mask(previous))&&pc(r)<=rem)total+=memo[row+1][r][rem-pc(r)];memo[row][previous][rem]=total;}
}
template<class F> static void enum_rec(int row,Mask previous,int rem,Mask full,F&& f){
  if(!suffix(row,previous,rem))return;
  if(row==rows){if(rem==0)f(full);return;}
  for(Mask r:valid)if(!(r&previous)&&pc(r)<=rem)enum_rec(row+1,r,rem-pc(r),full|(r<<((row+2)*n)),f);
}
template<class F> static void enum_layer(int k,F&& f){enum_rec(0,0,k,0,std::forward<F>(f));}
struct Block {int row,rem;Mask previous,full;u64 count,base;};
static void split_rec(int stop,int row,Mask previous,int rem,Mask full,std::vector<Block>& out){
  if(row==stop||row==rows){u64 z=suffix(row,previous,rem);if(z)out.push_back({row,rem,previous,full,z,0});return;}
  for(Mask r:valid)if(!(r&previous)&&pc(r)<=rem)split_rec(stop,row+1,r,rem-pc(r),full|(r<<((row+2)*n)),out);
}
static std::vector<Block> blocks_for(int k){std::vector<Block> b;split_rec(std::min(2,rows),0,0,k,0,b);u64 base=0;for(auto&x:b){x.base=base;base+=x.count;}if(base!=suffix(0,0,k))throw std::runtime_error("block census mismatch");return b;}
static u64 factorial_mod(int q){u64 a=1;for(int i=2;i<=q;i++)a=mulmod(a,u64(i));return a;}
static u64 payload_bytes(size_t previous,size_t current,u64 hash_bytes){u128 z=u128(previous+current)*8+hash_bytes;if(z>PAYLOAD_LIMIT_BYTES)throw std::runtime_error("payload guard: modeled layer allocation exceeds 7 GiB");return u64(z);}

static u64 compute_layer(int k,const std::vector<Block>& blocks,const std::vector<u64>& prev,const Hash& hprev,std::vector<u64>& cur,int threads){
  std::atomic<size_t> next{0};std::atomic<bool> bad{false};std::exception_ptr ep;std::mutex ep_mu;std::vector<u64> local(threads,0);std::vector<std::thread> pool;
  for(int t=0;t<threads;t++)pool.emplace_back([&,t]{try{while(!bad.load()){size_t q=next.fetch_add(1);if(q>=blocks.size())break;const Block& b=blocks[q];u64 at=b.base,subtotal=0;
    enum_rec(b.row,b.previous,b.rem,b.full,[&](Mask I){u128 sum=0;Mask r=I,cl=0;while(r){Mask bb=r&-r;r-=bb;cl|=closed[bit(bb)];sum+=(k==1?prev[0]:prev[hprev.find(I^bb)]);}int den=pc(cl);if(den<1||den>56)throw std::runtime_error("recurrence denominator outside field audit range");u64 value=mulmod(u64(sum%modp),inv[den]);cur[at++]=value;subtotal=addmod(subtotal,value);});
    if(at!=b.base+b.count)throw std::runtime_error("exclusive block write length mismatch");local[t]=addmod(local[t],subtotal);
  }}catch(...){bad.store(true);std::lock_guard<std::mutex> g(ep_mu);if(!ep)ep=std::current_exception();}});
  for(auto&th:pool)th.join();if(ep)std::rethrow_exception(ep);u64 total=0;for(u64 x:local)total=addmod(total,x);return total;
}

static void boundary_test(const std::string& out){
  Mask full=~Mask(0), top=(Mask(1)<<8)-1, high=Mask(1)<<63;std::array<Mask,64> a{};
  for(int v=0;v<64;v++){int x=v%8,y=v/8;for(auto [X,Y]:{std::pair{x-1,y},std::pair{x+1,y},std::pair{x,y-1},std::pair{x,y+1}})if(0<=X&&X<8&&0<=Y&&Y<8)a[v]|=Mask(1)<<(Y*8+X);a[v]|=Mask(1)<<v;}
  if((full&high)!=high||pc(top)!=8||pc(high)!=1||!(a[63]&high)||!(a[63]&(Mask(1)<<62))||!(a[63]&(Mask(1)<<55)))throw std::runtime_error("64-bit boundary failed");
  std::ofstream f(out);f<<"{\"N\":64,\"FULL_is_all_ones\":true,\"top_mask\":\"0xff\",\"bit63_nonzero\":true,\"corner63_closed_contains_55_62_63\":true}\n";
}
static void plan_n8(const std::string& out){
  constexpr size_t previous=159834860,current=162959262;u64 hash=Hash::bytes_for(previous),compute=payload_bytes(previous,current,hash),rebuild=payload_bytes(0,current,Hash::bytes_for(current));
  std::ofstream f(out);f<<"{\"n\":8,\"transition\":\"11-to-12\",\"previous\":"<<previous<<",\"current\":"<<current<<",\"hash_slots\":"<<Hash::slots_for(previous)<<",\"hash_load\":\"0.5954312533140182\",\"compute_payload_bytes\":"<<compute<<",\"rebuild_payload_bytes\":"<<rebuild<<",\"payload_limit_bytes\":"<<PAYLOAD_LIMIT_BYTES<<",\"guard_passes_without_allocating_layers\":true}\n";
}

int main(int argc,char**argv){try{
  if(argc==3&&std::string(argv[1])=="--boundary-test"){boundary_test(argv[2]);return 0;}
  if(argc==3&&std::string(argv[1])=="--plan-n8"){plan_n8(argv[2]);return 0;}
  if(argc<5||argc>6){std::cerr<<"usage: completion_layer_mod n prime_index threads output.json [--enable-n8]\n";return 2;}
  n=std::stoi(argv[1]);int pi=std::stoi(argv[2]),threads=std::stoi(argv[3]);std::string out=argv[4];bool enable8=argc==6&&std::string(argv[5])=="--enable-n8";
  if(n<3||n>8||pi<0||pi>=4||threads<1||threads>64||(n==8&&!enable8))return 2;if(std::filesystem::exists(out))throw std::runtime_error("refusing to overwrite residue checkpoint");
  modp=PRIMES[pi];for(int d=1;d<=64;d++)inv[d]=powmod(u64(d),modp-2);N=n*n;rows=n-2;alpha=(rows*n+1)/2;S=(Mask(1)<<n)-1;
  for(int v=0;v<N;v++){int x=v%n,y=v/n;for(auto [X,Y]:{std::pair{x-1,y},std::pair{x+1,y},std::pair{x,y-1},std::pair{x,y+1}})if(0<=X&&X<n&&0<=Y&&Y<n)closed[v]|=Mask(1)<<(Y*n+X);closed[v]|=Mask(1)<<v;}
  for(Mask r=0;r<(Mask(1)<<n);r++)if(!(r&(r<<1)))valid.push_back(r);std::memset(memo,0,sizeof(memo));precompute_suffix();
  int m=N-n;std::vector<u64> prev(1,factorial_mod(m));Hash hprev;u64 residue=prev[0],total_sets=1,max_layer=1,max_payload=8;
  for(int k=1;k<=alpha;k++){
    auto began=std::chrono::steady_clock::now();auto blocks=blocks_for(k);u64 amount=suffix(0,0,k);u64 oldcap=hprev.cap;u64 compute_payload=payload_bytes(prev.size(),size_t(amount),hprev.bytes());max_payload=std::max(max_payload,compute_payload);std::vector<u64> cur((size_t(amount)));u64 layer_sum=compute_layer(k,blocks,prev,hprev,cur,threads);
    if(k&1)residue=(residue+modp-layer_sum)%modp;else residue=addmod(residue,layer_sum);total_sets+=amount;max_layer=std::max(max_layer,amount);
    std::vector<u64>().swap(prev);hprev=Hash();u64 rebuild_payload=0,nextcap=0; if(k<alpha){rebuild_payload=payload_bytes(0,size_t(amount),Hash::bytes_for(size_t(amount)));max_payload=std::max(max_payload,rebuild_payload);Hash hcur{size_t(amount)};nextcap=hcur.cap;u64 at=0;enum_layer(k,[&](Mask I){hcur.insert(I,uint32_t(at++));});if(at!=amount)throw std::runtime_error("re-enumeration length mismatch");prev.swap(cur);hprev=std::move(hcur);}
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();std::cerr<<"layer="<<k<<" count="<<amount<<" old_hashcap="<<oldcap<<" compute_payload_bytes="<<compute_payload<<" next_hashcap="<<nextcap<<" rebuild_payload_bytes="<<rebuild_payload<<" seconds="<<seconds<<"\n";
  }
  std::string tmp=out+".tmp";{std::ofstream f(tmp);if(!f)throw std::runtime_error("cannot create checkpoint");f<<"{\"n\":"<<n<<",\"prime_index\":"<<pi<<",\"prime\":\""<<modp<<"\",\"threads\":"<<threads<<",\"algorithm\":\"modular-layer-rolling-reenumerated-hash\",\"residue\":\""<<residue<<"\",\"independent_sets\":"<<total_sets<<",\"max_layer_count\":"<<max_layer<<",\"payload_limit_bytes\":"<<PAYLOAD_LIMIT_BYTES<<",\"max_modeled_payload_bytes\":"<<max_payload<<"}\n";}std::filesystem::rename(tmp,out);return 0;
 }catch(const std::exception&e){std::cerr<<"completion_layer_mod: "<<e.what()<<"\n";return 1;}}
