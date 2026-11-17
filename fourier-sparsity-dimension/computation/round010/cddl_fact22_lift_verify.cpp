// Exact GF(2^14) / Walsh verifier for CDDL Fact 22 seed and triangular lifts.
// Coordinates are polynomial-basis bits x=sum x_i T^i for the emitted modulus.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using u32=uint32_t; using u64=uint64_t;
static int par(u32 x){return __builtin_parity(x);}
static u32 pmul(u32 a,u32 b){u32 r=0;while(b){if(b&1)r^=a;a<<=1;b>>=1;}return r;}
static int pdeg(u32 x){return x?31-__builtin_clz(x):-1;}
static u32 pmod(u32 a,u32 b){int d=pdeg(b);while(pdeg(a)>=d)a^=b<<(pdeg(a)-d);return a;}
static u32 pgcd(u32 a,u32 b){while(b){u32 t=pmod(a,b);a=b;b=t;}return a;}
struct GF {u32 mod; u32 mask=(1u<<14)-1;
 u32 mul(u32 a,u32 b)const{u32 r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&(1u<<14))a^=mod;}return r&mask;}
 u32 sq(u32 a)const{return mul(a,a);} u32 pow(u32 a,u32 e)const{u32 r=1;while(e){if(e&1)r=mul(r,a);a=mul(a,a);e>>=1;}return r;}
 u32 trace(u32 a)const{u32 r=a,t=a;for(int i=1;i<14;i++){t=sq(t);r^=t;} if(r!=0 && r!=1)throw std::runtime_error("trace not base field");return r;}
};
static bool irreducible(u32 p){ // Rabin test: x^(2^14)=x; q=2,7 prime divisors of 14.
 if(!(p&1) || pdeg(p)!=14)return false; u32 x=2, y=x;
 for(int i=0;i<14;i++)y=pmod(pmul(y,y),p); if(y!=x)return false;
 for(int e: {2,7}){y=x;for(int i=0;i<14/e;i++)y=pmod(pmul(y,y),p);if(pgcd(y^x,p)!=1)return false;}return true;
}
static u32 find_mod(){for(u32 low=1;low<(1u<<14);low+=2){u32 p=(1u<<14)|low;if(irreducible(p))return p;}throw std::runtime_error("no modulus");}
static void fwht(std::vector<int>& a){for(size_t h=1;h<a.size();h<<=1)for(size_t p=0;p<a.size();p+=2*h)for(size_t i=p;i<p+h;i++){int x=a[i],y=a[i+h];a[i]=x+y;a[i+h]=x-y;}}
static std::string hexmask(const std::vector<unsigned char>& f){std::ostringstream s;s<<"0x";for(int nib=(static_cast<int>(f.size())+3)/4-1;nib>=0;--nib){int v=0;for(int j=0;j<4;j++){int k=4*nib+j;if(k<static_cast<int>(f.size()))v|=f[k]<<j;}s<<"0123456789abcdef"[v];}return s.str();}
int main(int argc,char**argv){
 int maxl=4,samples=1;u64 seed=2026111710ULL;std::string out="cddl-lift.json";
 for(int i=1;i<argc;i++){std::string a=argv[i];auto take=[&](auto&v){if(++i>=argc)throw std::runtime_error("missing option");std::istringstream(argv[i])>>v;};if(a=="--max-l")take(maxl);else if(a=="--samples")take(samples);else if(a=="--seed")take(seed);else if(a=="--out")take(out);else if(a=="--help"){std::cout<<"--max-l L(0..4) --samples S --seed U64 --out FILE\n";return 0;}else throw std::runtime_error("unknown option");}
 if(maxl<0||maxl>4||samples<1)throw std::runtime_error("range"); auto begun=std::chrono::steady_clock::now();
 u32 mod=find_mod(); GF gf{mod};u32 alpha=0;for(u32 a=2;a<(1u<<14);a++)if((gf.sq(a)^a^1)==0){alpha=a;break;}if(!alpha)throw std::runtime_error("no F4 root");
 const int NX=1<<14;std::vector<unsigned char> f(NX);std::vector<int> sw(NX);for(int x=0;x<NX;x++){u32 xp=gf.pow((u32)x,57);f[x]=(unsigned char)gf.trace(gf.mul(alpha,xp));sw[x]=f[x]?-1:1;}auto Wseed=sw;fwht(Wseed);bool seedbent=true;int walsh_pos=0,walsh_neg=0;for(int z:Wseed){if(std::abs(z)!=128)seedbent=false;if(z==128)++walsh_pos;if(z==-128)++walsh_neg;}
 std::vector<unsigned char> anf=f;for(int b=0;b<14;b++)for(int x=0;x<NX;x++)if(x&(1<<b))anf[x]^=anf[x^(1<<b)];int deg=0;for(int x=0;x<NX;x++)if(anf[x])deg=std::max(deg,__builtin_popcount((u32)x));
 std::mt19937_64 rng(seed);std::ostringstream rec;rec<<"[";bool first=true,all_ok=seedbent;for(int l=0;l<=maxl;l++)for(int si=0;si<samples;si++){
  int Q=1<<l,n=14+2*l,N=1<<n,q=128*Q;std::vector<u32> lin(Q);std::vector<int> cc(Q),table(l);for(int w=0;w<Q;w++){lin[w]=rng()&((1u<<14)-1);cc[w]=rng()&1;}for(int j=0;j<l;j++)table[j]=rng()&((1<<(1<<j))-1);
  std::vector<int> pi(Q);bool bij=true;std::vector<int> seen(Q);for(int w=0;w<Q;w++){int aa=0;
   for(int j=0;j<l;j++)aa|=(((table[j]>>(w&((1<<j)-1)))&1)<<j);pi[w]=w^aa;if(seen[pi[w]]++)bij=false;}
  std::vector<int> vals(N);for(int idx=0;idx<N;idx++){int x=idx&((1<<14)-1),z=(idx>>14)&(Q-1),w=(idx>>(14+l))&(Q-1);int bit=f[x]^par((u32)(z&pi[w]))^par(lin[w]&(u32)x)^cc[w];vals[idx]=bit?-1:1;}
  auto W=vals;fwht(W);bool ok=bij;for(int v:W)if(std::abs(v)!=q){ok=false;break;}bool formula=true;std::vector<int> inv(Q);for(int w=0;w<Q;w++)inv[pi[w]]=w;for(int zeta=0;zeta<Q;zeta++)for(int omega=0;omega<Q;omega++){int w=inv[omega];for(int xi=0;xi<NX;xi++){int expected=Q*((par((u32)(zeta&w))^cc[w])?-1:1)*Wseed[xi^lin[w]];if(W[xi|(omega<<14)|(zeta<<(14+l))]!=expected)formula=false;}}
  auto bad=vals;bad[0]*=-1;fwht(bad);bool neg=false;for(int v:bad)if(std::abs(v)!=q){neg=true;break;}all_ok&=ok&&bij&&formula&&neg;
  if(!first)rec<<',';first=false;rec<<"{\"l\":"<<l<<",\"sample\":"<<si<<",\"variables\":"<<n<<",\"q\":"<<q<<",\"triangular_prefix_truth_tables\":[";for(int j=0;j<l;j++){if(j)rec<<',';rec<<table[j];}rec<<"],\"lin_by_w\":[";for(int w=0;w<Q;w++){if(w)rec<<',';rec<<lin[w];}rec<<"],\"constant_by_w\":[";for(int w=0;w<Q;w++){if(w)rec<<',';rec<<cc[w];}rec<<"],\"pi_bijective\":"<<(bij?"true":"false")<<",\"all_walsh_abs_q\":"<<(ok?"true":"false")<<",\"direct_walsh_matches_formula\":"<<(formula?"true":"false")<<",\"one_bit_flip_breaks_flatness\":"<<(neg?"true":"false")<<"}";
 }
 rec<<"]";double secs=std::chrono::duration<double>(std::chrono::steady_clock::now()-begun).count();std::ofstream os(out);
 os<<"{\n\"status\":\""<<(all_ok?"passed":"failed")<<"\",\"coordinate_convention\":\"polynomial-basis bits x_i*T^i, low bit is T^0\",\"degree14_modulus_hex\":\"0x"<<std::hex<<mod<<std::dec<<"\",\"rabin_irreducible\":true,\"alpha_hex\":\"0x"<<std::hex<<alpha<<std::dec<<"\",\"alpha_square_plus_alpha_plus_one\":true,\"seed_variables\":14,\"seed_walsh_all_abs_128\":"<<(seedbent?"true":"false")<<",\"seed_walsh_histogram\":{\"-128\":"<<walsh_neg<<",\"128\":"<<walsh_pos<<"},\"seed_anf_degree\":"<<deg<<",\"seed_truthmask_hex_le\":\""<<hexmask(f)<<"\",\"seed\":"<<seed<<",\"samples\":"<<rec.str()<<",\"elapsed_seconds\":"<<secs<<"\n}\n";
 std::cout<<"{\"status\":\""<<(all_ok?"passed":"failed")<<"\",\"modulus\":"<<mod<<",\"alpha\":"<<alpha<<",\"anf_degree\":"<<deg<<",\"elapsed_seconds\":"<<secs<<"}\n";return all_ok?0:2;
}
