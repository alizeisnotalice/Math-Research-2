// Exact joint (L, Manhattan distance, squared Euclidean distance) table.
// Standalone two-stage subset DP on a one-vertex-punctured square grid.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using u128=unsigned __int128; using u32=std::uint32_t;
static constexpr u128 UMAX=~u128(0);
static std::string dec(u128 x){if(!x)return "0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
static u128 addc(u128 a,u128 b){if(UMAX-a<b)throw std::overflow_error("u128 add");return a+b;}
static u128 mulc(u128 a,u128 b){if(a&&b>UMAX/a)throw std::overflow_error("u128 mul");return a*b;}
struct Grid{
 int n,N,M,hole; std::vector<int> orig,local; std::vector<u32> adj; u32 nv=0,full=0;
 Grid(int n_,int v):n(n_),N(n*n),M(N-1),hole(v),local(N,-1){
  if(n<2||M>31||v<0||v>=N)throw std::runtime_error("bad arguments");
  for(int q=0;q<N;++q)if(q!=hole){local[q]=int(orig.size());orig.push_back(q);} adj.assign(M,0);full=(u32(1)<<M)-1;
  const int dx[4]={1,-1,0,0},dy[4]={0,0,1,-1};
  for(int i=0;i<M;++i){int q=orig[i],x=q%n,y=q/n;for(int k=0;k<4;++k){int X=x+dx[k],Y=y+dy[k];if(X<0||X>=n||Y<0||Y>=n)continue;int r=Y*n+X;if(r==hole)nv|=u32(1)<<i;else adj[i]|=u32(1)<<local[r];}}
 }
};
struct Result{int n,v,M,dmax,e2max;std::vector<u128> h,joint,dist;u128 B=0,A=0,E2=0;};
static std::size_t ix(const Result&r,int t,int d,int e){return (std::size_t(t)*(r.dmax+1)+d)*(r.e2max+1)+e;}
static Result compute(int n,int v){
 Grid g(n,v); const std::size_t states=std::size_t(1)<<g.M; std::vector<u32> nb(states);for(u32 S=1;S<=g.full;++S){u32 b=S&(~S+1);nb[S]=nb[S^b]|g.adj[__builtin_ctz(b)];}
 std::vector<u128> f(states);f[g.full]=1;for(u32 S=g.full;S-->0;){u32 a=(g.nv|nb[S])&(g.full^S);u128 z=0;while(a){u32 b=a&(~a+1);a^=b;z=addc(z,f[S|b]);}f[S]=z;}
 Result r{};r.n=n;r.v=v;r.M=g.M;r.dmax=2*n-2;r.e2max=2*(n-1)*(n-1);r.h.assign(g.M+1,0);r.joint.assign(std::size_t(g.M+1)*(r.dmax+1)*(r.e2max+1),0);r.dist.assign(r.dmax+1,0);
 const u32 exterior=g.full&~g.nv;
 for(u32 qq=exterior;qq;){u32 bp=qq&(~qq+1);qq^=bp;int p=__builtin_ctz(bp);u32 allowed=g.full&~(bp|g.adj[p]);std::vector<int> verts;for(u32 a=allowed;a;){u32 b=a&(~a+1);a^=b;verts.push_back(__builtin_ctz(b));}
  int K=int(verts.size());std::vector<u128> ways(std::size_t(1)<<K);ways[0]=1;std::vector<u32> expand(std::size_t(1)<<K);expand[0]=bp;for(u32 I=1;I<(u32(1)<<K);++I){u32 b=I&(~I+1);expand[I]=expand[I^b]|(u32(1)<<verts[__builtin_ctz(b)]);}
  std::vector<u128> hp(g.M+1,0);for(u32 I=0;I<(u32(1)<<K);++I){if(!ways[I])continue;u32 S=expand[I];u128 suffix=0;for(u32 a=g.adj[p];a;){u32 b=a&(~a+1);a^=b;suffix=addc(suffix,f[S|b]);}int t=__builtin_popcount(S);hp[t]=addc(hp[t],mulc(ways[I],suffix));u32 elig=(g.nv|nb[S])&allowed&~S;while(elig){u32 b=elig&(~elig+1);elig^=b;int x=__builtin_ctz(b),j=0;while(verts[j]!=x)++j;ways[I|(u32(1)<<j)]=addc(ways[I|(u32(1)<<j)],ways[I]);}}
  u128 chk=0;for(int t=1;t<=g.M;++t)chk=addc(chk,hp[t]);if(chk!=f[bp])throw std::runtime_error("fibre partition");
  int po=g.orig[p],vx=v%n,vy=v/n,px=po%n,py=po/n,d=std::abs(px-vx)+std::abs(py-vy),e=(px-vx)*(px-vx)+(py-vy)*(py-vy);
  for(int t=1;t<=g.M;++t)if(hp[t]){r.h[t]=addc(r.h[t],hp[t]);r.joint[ix(r,t,d,e)]=addc(r.joint[ix(r,t,d,e)],hp[t]);}
 }
 for(int t=1;t<=g.M;++t)for(int d=0;d<=r.dmax;++d)for(int e=0;e<=r.e2max;++e){u128 z=r.joint[ix(r,t,d,e)];if(!z)continue;u128 w=mulc(t,z);r.B=addc(r.B,z);r.A=addc(r.A,w);r.dist[d]=addc(r.dist[d],w);r.E2=addc(r.E2,mulc(e,w));}
 return r;
}
static void emit(const Result&r,bool comma){std::cout<<(comma?",\n":"")<<"    {\"v\":"<<r.v<<",\"xy\":["<<r.v%r.n<<","<<r.v/r.n<<"],\"B_v_from_D\":\""<<dec(r.B)<<"\",\"A2_v_from_fibres\":\""<<dec(r.A)<<"\",\"weighted_e2_total\":\""<<dec(r.E2)<<"\",\"histogram_L\":[";for(int t=1;t<=r.M;++t){if(t>1)std::cout<<',';std::cout<<'"'<<dec(r.h[t])<<'"';}std::cout<<"],\"weighted_distance_counts\":[";for(int d=0;d<=r.dmax;++d){if(d)std::cout<<',';std::cout<<'"'<<dec(r.dist[d])<<'"';}std::cout<<"],\"joint_H_L_d_e2\":[";bool first=true;for(int t=1;t<=r.M;++t)for(int d=0;d<=r.dmax;++d)for(int e=0;e<=r.e2max;++e){u128 z=r.joint[ix(r,t,d,e)];if(!z)continue;if(!first)std::cout<<',';first=false;std::cout<<"["<<t<<","<<d<<","<<e<<",\""<<dec(z)<<"\"]";}std::cout<<"]}";}
int main(int argc,char**argv){if(argc<3){std::cerr<<"usage: gap_geometry_dp n v|all\n";return 2;}try{int n=std::stoi(argv[1]);std::vector<int> vs;if(std::string(argv[2])=="all")for(int v=0;v<n*n;++v)vs.push_back(v);else vs.push_back(std::stoi(argv[2]));std::cout<<"{\"n\":"<<n<<",\"method\":\"exact shared-f suffix plus compressed prefix DP\",\"u128_checked\":true,\"holes\":[";bool comma=false;for(int v:vs){emit(compute(n,v),comma);comma=true;}std::cout<<"]}\n";}catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<"\n";return 1;}}
