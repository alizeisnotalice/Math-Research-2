// Finite exact crosscheck for the R20 common-coordinate MM projection proof.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using U=std::uint64_t;

static int parity(int x){return __builtin_parity((unsigned)x);}
static int rank(std::vector<int> v){int p[16]={};int r=0;for(int x:v)while(x){int b=31-__builtin_clz((unsigned)x);if(p[b])x^=p[b];else{p[b]=x;++r;break;}}return r;}
static std::vector<int> span_basis(const std::vector<int>& candidates,int d){std::vector<int> out;for(int x:candidates)if(rank([&]{auto z=out;z.push_back(x);return z;}())>(int)out.size()){out.push_back(x);if((int)out.size()==d)break;}if((int)out.size()!=d)throw std::runtime_error("basis construction");return out;}
static std::vector<int> kernel_one(int lambda,int k){std::vector<int> c;for(int x=1;x<(1<<k);++x)if(!parity(x&lambda))c.push_back(x);return span_basis(c,k-1);}
static std::vector<int> kernel_two(int l0,int l1,int k){std::vector<int> c;for(int x=1;x<(1<<k);++x)if(!parity(x&l0)&&!parity(x&l1))c.push_back(x);return span_basis(c,k-2);}
static void fwht(std::vector<int>& a){for(int h=1;h<(int)a.size();h<<=1)for(int p=0;p<(int)a.size();p+=2*h)for(int i=p;i<p+h;++i){int u=a[i],v=a[i+h];a[i]=u+v;a[i+h]=u-v;}}
static bool get(const std::vector<U>& m,int i){return (m[i>>6]>>(i&63))&1ULL;}
static void set(std::vector<U>&m,int i){m[i>>6]|=1ULL<<(i&63);}
static std::string hx(U x){std::ostringstream o;o<<"0x"<<std::hex<<x;return o.str();}
static std::string words_json(const std::vector<U>& v){std::ostringstream o;o<<"[";for(size_t i=0;i<v.size();++i){if(i)o<<",";o<<"\""<<hx(v[i])<<"\"";}return o.str()+"]";}

struct Rec{int k,id;std::string kind;std::vector<U>f,g;std::vector<int>kbasis;int support=0;bool bentf=false,bentg=false,injective=false;};
static std::vector<U> truth(int k,const std::vector<int>& p,const std::vector<int>& h){int q=1<<k,N=q*q;std::vector<U> out((N+63)/64);for(int y=0;y<q;++y)for(int x=0;x<q;++x)if(parity(x&p[y])^h[y])set(out,x|(y<<k));return out;}
static bool bent(const std::vector<U>& f,int q){int N=q*q;std::vector<int>a(N);for(int i=0;i<N;++i)a[i]=get(f,i)?-1:1;fwht(a);for(int z:a)if(std::abs(z)!=q)return false;return true;}
static std::vector<U> xor_masks(const std::vector<U>&a,const std::vector<U>&b){auto c=a;for(size_t i=0;i<c.size();++i)c[i]^=b[i];return c;}
static bool injective_kernel(const std::vector<U>& A,const std::vector<int>&basis,int k){int N=1<<(2*k);for(int m=1;m<(1<<(int)basis.size());++m){int h=0;for(size_t i=0;i<basis.size();++i)if(m>>i&1)h^=basis[i];for(int z=0;z<N;++z)if(get(A,z)&&get(A,z^h))return false;}return true;}
static Rec make_case(int k,int id,const std::string&kind,std::mt19937_64& rng){int q=1<<k;std::vector<int> p(q),s(q),h(q),j(q),dom(q);for(int i=0;i<q;++i){p[i]=dom[i]=i;h[i]=int(rng()&1);}std::shuffle(p.begin(),p.end(),rng);std::shuffle(dom.begin(),dom.end(),rng);s=p;std::vector<int> r(q);std::vector<int> kb;
  if(kind=="c4_t0"){
    std::vector<int>S(dom.begin(),dom.begin()+4),perm{0,1,2,3};do{std::shuffle(perm.begin(),perm.end(),rng);}while(std::any_of(perm.begin(),perm.end(),[i=0](int z)mutable{return z==i++;}));
    for(int i=0;i<4;++i){s[S[i]]=p[S[perm[i]]];r[S[i]]=int(rng()&1);}
    std::vector<int>d;for(int y:S)d.push_back(p[y]^s[y]);int v=-1;for(int x=0;x<q;++x){bool ok=true;for(int z:d)if(!parity(z&x)){ok=false;break;}if(ok){v=x;break;}}if(v<0)throw std::runtime_error("c4 v");
    int l0=-1,l1=-1;for(int a=1;a<q&&l0<0;++a)for(int b=a+1;b<q;++b)if(rank({a,b})==2){bool distinct=true;std::set<int> im;for(int y:S)im.insert(parity(a&y)|(parity(b&y)<<1));if(im.size()==4){l0=a;l1=b;break;}}if(l0<0)throw std::runtime_error("c4 quotient");
    for(int w:kernel_two(l0,l1,k))kb.push_back(w<<k);kb.push_back(v);
  }else if(kind=="c2_t1"){
    int a=dom[0],b=dom[1],c=dom[2];s[a]=p[b];s[b]=p[a];r[a]=int(rng()&1);r[b]=int(rng()&1);r[c]=1;int d=p[a]^p[b],ea=1^r[a],eb=1^r[b],lam=-1;for(int z=1;z<q;++z)if(parity(z&(c^a))&&parity(z&(c^b))){lam=z;break;}if(lam<0)throw std::runtime_error("c2 lambda");
    int delta=a^b;auto W=kernel_one(lam,k);std::vector<int> candidates{delta};for(int x:W)candidates.push_back(x);W=span_basis(candidates,k-1);int rhs=1^ea^eb,u=0;for(int x=0;x<q;++x)if(parity(d&x)==rhs){u=x;break;}for(size_t i=0;i<W.size();++i)kb.push_back(((i==0?u:0))|(W[i]<<k));
  }else if(kind=="c0_t2"){
    int a=dom[0],b=dom[1];r[a]=r[b]=1;int delta=a^b,lam=0;for(int z=1;z<q;++z)if(parity(z&delta)){lam=z;break;}for(int w:kernel_one(lam,k))kb.push_back(w<<k);
  }else throw std::runtime_error("kind");
  for(int y=0;y<q;++y)j[y]=h[y]^r[y];Rec out{k,id,kind,truth(k,p,h),truth(k,s,j),kb};auto A=xor_masks(out.f,out.g);for(U w:A)out.support+=__builtin_popcountll(w);out.bentf=bent(out.f,q);out.bentg=bent(out.g,q);out.injective=injective_kernel(A,kb,k);if(out.support!=2*q||rank(kb)!=k-1||!out.bentf||!out.bentg||!out.injective)throw std::runtime_error("case verification");return out;}

int main(int ac,char**av){if(ac!=3){std::cerr<<"usage: verify_common_mm_pairs smoke|full output.json\n";return 2;}std::string mode=av[1];if(mode!="smoke"&&mode!="full")return 2;auto began=std::chrono::steady_clock::now();std::mt19937_64 rng(20261120ULL);std::vector<Rec> rows;std::array<std::string,3> kinds{"c4_t0","c2_t1","c0_t2"};int per=(mode=="smoke"?1:8);for(int k=2;k<=6;++k)for(const auto&kind:kinds)for(int id=0;id<per;++id)rows.push_back(make_case(k,id,kind,rng));
  std::ofstream f(av[2]);if(!f)return 3;f<<"{\n\"schema\":\"a28-r020-common-mm-crosscheck-v1\",\n\"mode\":\""<<mode<<"\",\n\"seed\":20261120,\n\"input_indexing\":\"index=x+(y<<k), x and y each k-bit, little-endian 64-bit truth words\",\n\"cases\":[";for(size_t i=0;i<rows.size();++i){auto&r=rows[i];if(i)f<<",";f<<"{\"k\":"<<r.k<<",\"n\":"<<2*r.k<<",\"N\":"<<(1<<(2*r.k))<<",\"type\":\""<<r.kind<<"\",\"sample\":"<<r.id<<",\"support_size\":"<<r.support<<",\"kernel_basis_ambient\":[";for(size_t z=0;z<r.kbasis.size();++z){if(z)f<<",";f<<r.kbasis[z];}f<<"],\"f_bent\":true,\"g_bent\":true,\"projection_injective\":true,\"f_truth_words_le\":"<<words_json(r.f)<<",\"g_truth_words_le\":"<<words_json(r.g)<<"}";}f<<"],\n\"case_count\":"<<rows.size()<<",\n\"all_checks_passed\":true,\n\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<"\n}\n";return 0;}
