// Exact external-memory inclusion--exclusion calculation for the L-boundary
// family.  The program deliberately does not share code with the earlier n=8
// verifier: it enumerates row states, quotients only by diagonal reflection,
// and reconstructs the full alternating sums orbit by orbit.
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace fs = std::filesystem;
constexpr int CW = 5;                    // 320 bit C(I)
constexpr int AW = 6;                    // 384 bit alternating accumulators
using U64 = uint64_t;

struct U320 { std::array<U64,CW> w{}; };
struct U384 { std::array<U64,AW> w{}; };

static void add320(U320& a, const U320& b) {
  unsigned __int128 c=0;
  for (int i=0;i<CW;i++) { c=(unsigned __int128)a.w[i]+b.w[i]+(U64)c; a.w[i]=(U64)c; c>>=64; }
  assert(c==0);
}
static void div320(U320& a, unsigned d) {
  unsigned __int128 r=0;
  for (int i=CW-1;i>=0;i--) { unsigned __int128 x=(r<<64)|a.w[i]; a.w[i]=(U64)(x/d); r=x%d; }
  assert(r==0);
}
static void add384(U384& a, const U384& b) {
  unsigned __int128 c=0;
  for (int i=0;i<AW;i++) { c=(unsigned __int128)a.w[i]+b.w[i]+(U64)c; a.w[i]=(U64)c; c>>=64; }
  assert(c==0);
}
static void sub384(U384& a, const U384& b) {
  U64 borrow=0;
  for (int i=0;i<AW;i++) { U64 x=b.w[i]+borrow; U64 next=(x<b.w[i]); U64 old=a.w[i]; a.w[i]=old-x; borrow=(old<x)|next; }
  assert(!borrow);
}
static U384 mulsmall(const U320& x, unsigned a) {
  U384 z; unsigned __int128 c=0;
  for (int i=0;i<CW;i++) { c=(unsigned __int128)x.w[i]*a+(U64)c; z.w[i]=(U64)c; c>>=64; }
  z.w[CW]=(U64)c; return z;
}
static U384 u320to384(const U320& x) { U384 z; for(int i=0;i<CW;i++)z.w[i]=x.w[i]; return z; }
static bool divisible(const U384& x, U64 d) { U64 r=0; for(int i=AW-1;i>=0;i--){ unsigned __int128 q=((unsigned __int128)r<<64)|x.w[i]; r=(U64)(q%d); } return r==0; }
static U384 div384(U384 x, U64 d) { U64 r=0; for(int i=AW-1;i>=0;i--){ unsigned __int128 q=((unsigned __int128)r<<64)|x.w[i]; x.w[i]=(U64)(q/d); r=(U64)(q%d); } assert(r==0); return x; }
static std::string dec(U384 x) {
  bool zero=true; for(U64 v:x.w) zero&=(v==0); if(zero)return "0";
  std::string s;
  while(true) { U64 r=0; bool nz=false; for(int i=AW-1;i>=0;i--){ unsigned __int128 q=((unsigned __int128)r<<64)|x.w[i]; x.w[i]=(U64)(q/10); r=(U64)(q%10); nz|=(x.w[i]!=0); } s.push_back(char('0'+r)); if(!nz)break; }
  std::reverse(s.begin(),s.end()); return s;
}
static U320 factorial(int n) { U320 z; z.w[0]=1; for(int a=2;a<=n;a++){ unsigned __int128 c=0; for(int i=0;i<CW;i++){c=(unsigned __int128)z.w[i]*a+(U64)c;z.w[i]=(U64)c;c>>=64;} assert(!c); } return z; }
static U64 mix(U64 x) { x^=x>>30; x*=0xbf58476d1ce4e5b9ULL; x^=x>>27; x*=0x94d049bb133111ebULL; return x^(x>>31); }

struct Geometry {
  int n,q,m=0; U64 Ubits=0;
  bool S[9][9]{}, R[9][9]{}; int rid[9][9]{};
  U64 closedR[49]{}, rbit[49]{}; int reflpos[49]{};
  std::vector<U64> pairForbid;
  Geometry(int nn):n(nn),q(nn-2) {
    assert(n>=4 && n<=9 && q<=7);
    for(int x=2;x<n;x++) S[x][0]=true;
    for(int x=0;x<=2;x++) S[x][1]=true;
    for(int y=2;y<n;y++) S[0][y]=true;
    bool closedS[9][9]{};
    const int dx[5]={0,1,-1,0,0}, dy[5]={0,0,0,1,-1};
    for(int y=0;y<n;y++)for(int x=0;x<n;x++) if(S[x][y]) for(int a=0;a<5;a++){int X=x+dx[a],Y=y+dy[a];if(X>=0&&X<n&&Y>=0&&Y<n)closedS[X][Y]=true;}
    for(int y=0;y<n;y++)for(int x=0;x<n;x++)if(!S[x][y]) { R[x][y]=true; rid[x][y]=m++; }
    assert(m==n*n-(2*n-1));
    for(int y=0;y<q;y++)for(int x=0;x<q;x++) if(!(x==0&&y==0)) Ubits|=U64(1)<<(y*q+x);
    for(int y=0;y<q;y++)for(int x=0;x<q;x++) {
      int p=y*q+x; reflpos[p]=x*q+y;
      if(x==0&&y==0)continue;
      int X=x+2,Y=y+2; assert(!closedS[X][Y]);
      U64 c=0;
      for(int a=0;a<5;a++){int XX=X+dx[a],YY=Y+dy[a];if(XX<0||XX>=n||YY<0||YY>=n)continue;assert(R[XX][YY]);c|=U64(1)<<rid[XX][YY];}
      closedR[p]=c; rbit[p]=U64(1)<<p;
    }
    std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> pairs;
    bool boundary[9][9]{};
    for(int y=0;y<n;y++)for(int x=0;x<n;x++)if(S[x][y])for(int a=1;a<5;a++){int X=x+dx[a],Y=y+dy[a];if(X>=0&&X<n&&Y>=0&&Y<n&&!S[X][Y])boundary[X][Y]=true;}
    int b=0;
    for(int y=0;y<n;y++)for(int x=0;x<n;x++)if(boundary[x][y]){ b++; for(int a=1;a<5;a++){int X=x+dx[a],Y=y+dy[a];if(X>=0&&X<n&&Y>=0&&Y<n&&!closedS[X][Y]) pairs.emplace_back(std::make_pair(x,y),std::make_pair(X,Y));} }
    assert(b==2*n-2);
    for(auto [w,y]:pairs){ U64 f=0; for(auto z:{w,y}) for(int a=0;a<5;a++){int X=z.first+dx[a],Y=z.second+dy[a]; if(X>=2&&X<n&&Y>=2&&Y<n && !(X==2&&Y==2)) f|=U64(1)<<((Y-2)*q+(X-2));} pairForbid.push_back(f); }
  }
  U64 reflect(U64 I) const { U64 J=0; while(I){int p=__builtin_ctzll(I);I&=I-1;J|=U64(1)<<reflpos[p];}return J; }
  U64 canon(U64 I) const { U64 J=reflect(I); return std::min(I,J); }
  unsigned closedSize(U64 I) const { U64 r=0;while(I){int p=__builtin_ctzll(I);I&=I-1;r|=closedR[p];}return __builtin_popcountll(r); }
  unsigned eligible(U64 I) const { unsigned a=0;for(U64 f:pairForbid)a+=(I&f)==0;return a; }
};

static std::string kfile(const fs::path& root, int k, const char* ext="u64") { char b[80]; std::snprintf(b,sizeof b,"layers/k%02d.%s",k,ext); return (root/b).string(); }
static void writeall(int fd,const void* p,size_t n) { const char* q=(const char*)p;while(n){ssize_t z=::write(fd,q,n);assert(z>0);q+=z;n-=z;} }
static void pwriteall(int fd,const void* p,size_t n,off_t at) { const char* q=(const char*)p;while(n){ssize_t z=::pwrite(fd,q,n,at);assert(z>0);q+=z;n-=z;at+=z;} }

struct Generator {
  const Geometry& g; fs::path root; int nth; std::vector<U64> rows, first;
  Generator(const Geometry& gg,fs::path rr,int nt):g(gg),root(rr),nth(nt) { for(U64 r=0;r<(U64(1)<<g.q);r++)if(!(r&(r<<1)))rows.push_back(r);for(U64 r:rows)if(!(r&1))first.push_back(r); }
  void rec(int tid,int y,U64 mask,U64 prev,std::vector<std::vector<U64>>& buf,std::vector<int>& fd) {
    if(y==g.q){U64 J=g.reflect(mask);if(mask<=J){int k=__builtin_popcountll(mask);buf[k].push_back(mask);if(buf[k].size()>=16384){writeall(fd[k],buf[k].data(),buf[k].size()*8);buf[k].clear();}}return;}
    for(U64 r:rows)if(!(r&prev))rec(tid,y+1,mask|(r<<(y*g.q)),r,buf,fd);
  }
  void run() {
    fs::create_directories(root/"layers"); fs::create_directories(root/"segments");
    std::vector<std::pair<U64,U64>> seeds;for(U64 a:first)for(U64 b:rows)if(!(a&b))seeds.push_back({a,b});
    std::vector<std::thread> th;
    for(int t=0;t<nth;t++)th.emplace_back([&,t]{std::vector<std::vector<U64>> buf(g.q*g.q+1);std::vector<int> fd(g.q*g.q+1);for(int k=0;k<=g.q*g.q;k++){auto p=root/"segments"/("t"+std::to_string(t)+"-k"+(k<10?"0":"")+std::to_string(k)+".u64");fd[k]=::open(p.c_str(),O_CREAT|O_TRUNC|O_WRONLY,0644);assert(fd[k]>=0);}for(size_t j=t;j<seeds.size();j+=nth)rec(t,2,seeds[j].first|(seeds[j].second<<g.q),seeds[j].second,buf,fd);for(int k=0;k<=g.q*g.q;k++){if(!buf[k].empty())writeall(fd[k],buf[k].data(),buf[k].size()*8);::close(fd[k]);}});
    for(auto& x:th)x.join();
    // The empty set is absent from two-row prefixes.
    {int fd=::open(kfile(root,0).c_str(),O_CREAT|O_TRUNC|O_WRONLY,0644);U64 z=0;writeall(fd,&z,8);::close(fd);}
    for(int k=1;k<=g.q*g.q;k++){int out=::open(kfile(root,k).c_str(),O_CREAT|O_TRUNC|O_WRONLY,0644);assert(out>=0);std::vector<char>b(1<<20);for(int t=0;t<nth;t++){auto p=root/"segments"/("t"+std::to_string(t)+"-k"+(k<10?"0":"")+std::to_string(k)+".u64");int in=::open(p.c_str(),O_RDONLY);assert(in>=0);while(true){ssize_t z=::read(in,b.data(),b.size());assert(z>=0);if(!z)break;writeall(out,b.data(),z);}::close(in);}::close(out);}
  }
};

struct Hash {
  const U64* keys; size_t n,cap,mask; std::atomic<uint32_t>* slot;
  Hash(const U64* a,size_t z,int nth):keys(a),n(z){
    // Slots store index+1 in a uint32_t; zero remains the empty sentinel.
    // Do not use weak CAS here: a spurious failure can leave e==0, and the
    // old code then evaluated keys[e-1].
    assert(a && z && z<=size_t(UINT32_MAX-1) && nth>0);
    cap=1;
    while(cap<2*n) { assert(cap<=std::numeric_limits<size_t>::max()/2); cap<<=1; }
    mask=cap-1;
    slot=new std::atomic<uint32_t>[cap];
    std::vector<std::thread> th;
    for(int t=0;t<nth;t++)th.emplace_back([&,t]{for(size_t i=t;i<cap;i+=nth)slot[i].store(0,std::memory_order_relaxed);});
    for(auto&x:th)x.join(); th.clear();
    for(int t=0;t<nth;t++)th.emplace_back([&,t]{for(size_t i=t;i<n;i+=nth){
      size_t p=mix(keys[i])&mask;
      for(;;){
        uint32_t e=0;
        if(slot[p].compare_exchange_strong(e,uint32_t(i+1),std::memory_order_relaxed))break;
        if(e!=0 && keys[e-1]==keys[i]) { assert(e==i+1); break; }
        p=(p+1)&mask;
      }
    }});
    for(auto&x:th)x.join();
  }
  ~Hash(){delete[]slot;}
  size_t get(U64 key)const{size_t p=mix(key)&mask;for(;;){uint32_t e=slot[p].load(std::memory_order_relaxed);assert(e);if(keys[e-1]==key)return e-1;p=(p+1)&mask;}}
};

static size_t fsize(const std::string&p){struct stat s{};assert(!::stat(p.c_str(),&s));return (size_t)s.st_size;}
static void* mapread(const std::string&p,size_t n){int fd=::open(p.c_str(),O_RDONLY);assert(fd>=0);void*z=::mmap(nullptr,n,PROT_READ,MAP_SHARED,fd,0);assert(z!=MAP_FAILED);::close(fd);return z;}

struct Calculator {
  const Geometry& g; fs::path root; int nth; U384 hp{},hn{},tp{},tn{};
  Calculator(const Geometry&gg,fs::path rr,int nt):g(gg),root(rr),nth(nt){}
  void addterm(U64 I,const U320& c,unsigned orb,U384& lh,U384&lt) {
    U384 a=mulsmall(c,orb); U384 b=mulsmall(c,g.eligible(I)); if(orb==2)b=mulsmall(c,g.eligible(I)+g.eligible(g.reflect(I)));
    add384(lh,a);add384(lt,b);
  }
  void run() {
    U320 fact=factorial(g.m); hp=u320to384(fact); tp=mulsmall(fact,(unsigned)g.pairForbid.size());
    {int fd=::open(kfile(root,0,"u320").c_str(),O_CREAT|O_TRUNC|O_WRONLY,0644);writeall(fd,&fact,sizeof fact);::close(fd);}
    std::ofstream log(root/"calc-progress.log");log<<"m="<<g.m<<" pairs="<<g.pairForbid.size()<<"\n";
    for(int k=1;k<=g.q*g.q;k++) {
      auto tm=std::chrono::steady_clock::now();
      std::string pm=kfile(root,k-1), pv=kfile(root,k-1,"u320"), cm=kfile(root,k), cv=kfile(root,k,"u320");
      size_t pn=fsize(pm)/8, cn=fsize(cm)/8; if(cn==0) break; void* mm=mapread(pm,pn*8);void* vm=mapread(pv,pn*sizeof(U320));auto* masks=(const U64*)mm;auto* vals=(const U320*)vm;
      Hash hash(masks,pn,nth);
      int out=::open(cv.c_str(),O_CREAT|O_TRUNC|O_WRONLY,0644);assert(out>=0);assert(!::ftruncate(out,off_t(cn*sizeof(U320))));
      int infd=::open(cm.c_str(),O_RDONLY);assert(infd>=0);constexpr size_t B=32768;U384 kh{},kt{};
      for(size_t base=0;base<cn;base+=B){size_t z=std::min(B,cn-base);std::vector<U64> in(z);assert(::pread(infd,in.data(),z*8,base*8)==ssize_t(z*8));std::vector<U320> outv(z);std::vector<U384> hs(nth),ts(nth);std::vector<std::thread> th;
        for(int ti=0;ti<nth;ti++)th.emplace_back([&,ti]{for(size_t a=ti;a<z;a+=nth){U64 I=in[a], J=g.reflect(I);U320 sum{};for(U64 x=I;x;){int p=__builtin_ctzll(x);x&=x-1;U64 child=std::min(I^(U64(1)<<p),J^(U64(1)<<g.reflpos[p]));add320(sum,vals[hash.get(child)]);}div320(sum,g.closedSize(I));outv[a]=sum;addterm(I,sum,(I==J)?1:2,hs[ti],ts[ti]);}});
        for(auto&x:th)x.join();for(int ti=0;ti<nth;ti++){add384(kh,hs[ti]);add384(kt,ts[ti]);}pwriteall(out,outv.data(),z*sizeof(U320),off_t(base*sizeof(U320)));}
      ::close(infd);::close(out);::munmap(mm,pn*8);::munmap(vm,pn*sizeof(U320));assert(!::unlink(pv.c_str()));if(k&1){add384(hn,kh);add384(tn,kt);}else{add384(hp,kh);add384(tp,kt);}
      auto sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-tm).count();log<<"k="<<k<<" prev="<<pn<<" cur="<<cn<<" seconds="<<sec<<" hplus="<<dec(hp)<<" hminus="<<dec(hn)<<"\n";log.flush();
      std::cerr<<"calc k="<<k<<" prev="<<pn<<" cur="<<cn<<" sec="<<sec<<"\n";
    }
    U384 h=hp,t=tp;sub384(h,hn);sub384(t,tn);assert(divisible(t,U64(g.m)*(g.m-1)));U384 H=div384(t,U64(g.m)*(g.m-1));
    std::ofstream o(root/"result.json");o<<"{\n  \"n\": "<<g.n<<",\n  \"S\": "<<(2*g.n-1)<<",\n  \"m\": "<<g.m<<",\n  \"boundary\": "<<(2*g.n-2)<<",\n  \"U\": "<<(g.q*g.q-1)<<",\n  \"pairs\": "<<g.pairForbid.size()<<",\n  \"h\": \""<<dec(h)<<"\",\n  \"pair_numerator_over_m_m1\": \""<<dec(H)<<"\",\n  \"T\": \""<<dec(t)<<"\",\n  \"p_numerator\": \""<<dec(H)<<"\",\n  \"p_denominator\": \""<<dec(h)<<"\"\n}\n";
    std::cout<<"h="<<dec(h)<<"\nH="<<dec(H)<<"\nT="<<dec(t)<<"\n";
  }
};

static void audit_layers(const Geometry& g, const fs::path& root) {
  U64 raw_total=0, rep_total=0, fixed_total=0, maxrep=0; int maxk=-1;
  std::ofstream o(root/"reflection-layer-audit-cpp.json");
  o<<"{\n  \"n\": "<<g.n<<",\n  \"layers\": [\n";
  std::vector<U64> buf(1<<17);
  for(int k=0;k<=g.q*g.q;k++){
    std::string p=kfile(root,k); int fd=::open(p.c_str(),O_RDONLY);assert(fd>=0); U64 reps=0,fixed=0;
    while(true){ssize_t z=::read(fd,buf.data(),buf.size()*8);assert(z>=0);if(!z)break;assert(z%8==0);for(size_t i=0;i<size_t(z)/8;i++){U64 x=buf[i],rx=g.reflect(x);assert(x<=rx);reps++;fixed+=(x==rx);}}
    ::close(fd); U64 raw=2*reps-fixed;raw_total+=raw;rep_total+=reps;fixed_total+=fixed;if(reps>maxrep){maxrep=reps;maxk=k;}
    o<<"    {\"k\": "<<k<<", \"reps\": "<<reps<<", \"fixed\": "<<fixed<<", \"raw\": "<<raw<<"}"<<(k==g.q*g.q?"\n":" ,\n");
  }
  o<<"  ],\n  \"raw_total\": "<<raw_total<<",\n  \"rep_total\": "<<rep_total<<",\n  \"fixed_total\": "<<fixed_total<<",\n  \"max_reps\": "<<maxrep<<",\n  \"max_k\": "<<maxk<<"\n}\n";
  std::cout<<"audit raw_total="<<raw_total<<" rep_total="<<rep_total<<" fixed_total="<<fixed_total<<" max_reps="<<maxrep<<" max_k="<<maxk<<"\n";
}

int main(int argc,char**argv){
  if(argc!=5){std::cerr<<"usage: boundary_kernel N threads root generate|calculate|audit\n";return 2;}
  int n=std::stoi(argv[1]),nt=std::stoi(argv[2]);
  if(nt<=0){std::cerr<<"threads must be positive\n";return 2;}
  fs::path root=argv[3];std::string mode=argv[4];Geometry g(n);
  std::cerr<<"n="<<n<<" S="<<2*n-1<<" m="<<g.m<<" U="<<(g.q*g.q-1)<<" boundary="<<2*n-2<<" pairs="<<g.pairForbid.size()<<"\n";
  if(mode=="generate")Generator(g,root,nt).run();else if(mode=="calculate")Calculator(g,root,nt).run();else if(mode=="audit")audit_layers(g,root);else return 2;
}
