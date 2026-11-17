// Selected-prefix exact completion counts via independent-set inclusion-exclusion.
// Usage: completion_ie_selected n hex_mask [--cross]
// n<=7; all stored quantities use U256.  No floating arithmetic.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
using Mask=uint64_t;using u64=uint64_t;
struct U{std::array<u64,4>w{};};
static bool z(const U&a){return !(a.w[0]|a.w[1]|a.w[2]|a.w[3]);}
static int cmp(const U&a,const U&b){for(int i=3;i>=0;i--)if(a.w[i]!=b.w[i])return a.w[i]<b.w[i]?-1:1;return 0;}
static void add(U&a,const U&b){u64 c=0;for(int i=0;i<4;i++){u64 x=a.w[i]+b.w[i];bool c1=x<a.w[i];u64 y=x+c;bool c2=y<x;a.w[i]=y;c=c1||c2;}if(c)throw std::overflow_error("U256 add");}
static U sub(U a,const U&b){u64 br=0;for(int i=0;i<4;i++){u64 q=b.w[i]+br;bool c=q<b.w[i];u64 old=a.w[i];a.w[i]-=q;br=c||(old<q);}if(br)throw std::runtime_error("negative");return a;}
static U mul(U a,u64 q){U r;__uint128_t c=0;for(int i=0;i<4;i++){__uint128_t x=(__uint128_t)a.w[i]*q+c;r.w[i]=x;c=x>>64;}if(c)throw std::overflow_error("U256 mul");return r;}
static U div(U a,u64 q){U r;__uint128_t rem=0;for(int i=3;i>=0;i--){__uint128_t x=(rem<<64)|a.w[i];r.w[i]=x/q;rem=x%q;}if(rem)throw std::runtime_error("nonintegral");return r;}
static std::string dec(U a){if(z(a))return"0";std::string s;while(!z(a)){__uint128_t r=0;for(int i=3;i>=0;i--){__uint128_t x=(r<<64)|a.w[i];a.w[i]=x/10;r=x%10;}s.push_back('0'+int(r));}std::reverse(s.begin(),s.end());return s;}
static int pop(Mask x){return __builtin_popcountll(x);}static int bit(Mask x){return __builtin_ctzll(x);}
static int n,N;static Mask FULL;static std::array<Mask,49> adj,closed;
static Mask neigh(Mask s){Mask r=0;while(s){Mask b=s&-s;s-=b;r|=adj[bit(b)];}return r;}
static U factorial(int m){U a;a.w[0]=1;for(int i=2;i<=m;i++)a=mul(a,i);return a;}
static void enum_ind(Mask avail,Mask chosen,std::vector<Mask>&v){if(!avail){v.push_back(chosen);return;}Mask b=avail&-avail;int x=bit(b);enum_ind(avail^b,chosen,v);enum_ind((avail^b)&~adj[x],chosen|b,v);}
static std::string cells(Mask s){std::string o="[";bool f=1;for(int v=0;v<N;v++)if(s&(Mask(1)<<v)){if(!f)o+=",";f=0;o+="["+std::to_string(v%n)+","+std::to_string(v/n)+"]";}return o+"]";}

struct Store{std::vector<Mask> m;std::vector<U> c;std::unordered_map<Mask,uint32_t> where;};
static Store build_C(Mask S){
 Mask Uset=FULL&~neigh(S)&~S;Store q;enum_ind(Uset,0,q.m);std::sort(q.m.begin(),q.m.end(),[](Mask a,Mask b){return pop(a)<pop(b)||(pop(a)==pop(b)&&a<b);});q.c.resize(q.m.size());q.where.reserve(q.m.size()*2);for(uint32_t i=0;i<q.m.size();i++)q.where.emplace(q.m[i],i);
 int m=N-pop(S);q.c[0]=factorial(m);std::vector<size_t> begin(N+2,q.m.size());size_t at=0;for(int k=0;k<=N;k++){begin[k]=at;while(at<q.m.size()&&pop(q.m[at])==k)at++;}begin[N+1]=q.m.size();
 for(int k=1;k<=N;k++){size_t lo=begin[k],hi=(k==N?q.m.size():begin[k+1]);if(lo==hi)continue;std::vector<std::thread> th;for(int t=0;t<12;t++)th.emplace_back([&,t]{size_t a=lo+(hi-lo)*size_t(t)/12,b=lo+(hi-lo)*size_t(t+1)/12;for(size_t j=a;j<b;j++){Mask I=q.m[j],r=I,cl=0;U sum;while(r){Mask bb=r&-r;r-=bb;cl|=closed[bit(bb)];add(sum,q.c[q.where.at(I^bb)]);}q.c[j]=div(sum,pop(cl));}});for(auto&t:th)t.join();}
 return q;
}
static U F(const Store&q,Mask allowed){U pos,neg;for(size_t i=0;i<q.m.size();i++)if(!(q.m[i]&~allowed)){if(pop(q.m[i])&1)add(neg,q.c[i]);else add(pos,q.c[i]);}return sub(pos,neg);}

// Cross-check by the defining completion recurrence, over subsets of V\S.
static U direct_h(Mask S){std::vector<int> out;for(int v=0;v<N;v++)if(!(S&(Mask(1)<<v)))out.push_back(v);int m=out.size();if(m>23)throw std::runtime_error("cross complement too large");std::vector<int> ix(N,-1);for(int i=0;i<m;i++)ix[out[i]]=i;size_t M=size_t(1)<<m;std::vector<U> H(M);H[M-1].w[0]=1;for(size_t t=M-1;t-->0;){Mask cur=S;for(int i=0;i<m;i++)if(t&(size_t(1)<<i))cur|=Mask(1)<<out[i];Mask bd=neigh(cur)&~cur;U a;while(bd){Mask b=bd&-bd;bd-=b;int k=ix[bit(b)];add(a,H[t|(size_t(1)<<k)]);}H[t]=a;}return H[0];}

int main(int argc,char**argv){try{
 if(argc<3||argc>4)return 2;n=std::atoi(argv[1]);if(n<2||n>7)return 2;N=n*n;FULL=(Mask(1)<<N)-1;Mask S=std::stoull(argv[2],nullptr,0);if(S&~FULL||!S)return 2;bool cross=argc==4&&std::string(argv[3])=="--cross";
 for(int v=0;v<N;v++){int x=v%n,y=v/n;for(auto [X,Y]:{std::pair{x-1,y},std::pair{x+1,y},std::pair{x,y-1},std::pair{x,y+1}})if(0<=X&&X<n&&0<=Y&&Y<n)adj[v]|=Mask(1)<<(Y*n+X);closed[v]=adj[v]|(Mask(1)<<v);}int m=N-pop(S);if(m<2)return 2;Store q=build_C(S);Mask Uset=FULL&~neigh(S)&~S;U hs=F(q,Uset),num;Mask bd=neigh(S)&~S;int B=pop(bd),pairs=0;
 for(Mask a=bd;a;){Mask xb=a&-a;a-=xb;Mask ys=adj[bit(xb)]&~(S|neigh(S));while(ys){Mask yb=ys&-ys;ys-=yb;pairs++;Mask allow=FULL&~neigh(S|xb|yb)&~(S|xb|yb);U f=F(q,allow);add(num,div(f,u64(m)*(m-1)));}}
 U lhs=mul(num,B),rhs=mul(hs,4);std::cout<<"{\"n\":"<<n<<",\"S_cells\":"<<cells(S)<<",\"m\":"<<m<<",\"independent_sets\":"<<q.m.size()<<",\"boundary\":"<<B<<",\"dependent_pairs\":"<<pairs<<",\"hS\":\""<<dec(hs)<<"\",\"p_numerator\":\""<<dec(num)<<"\",\"p_denominator\":\""<<dec(hs)<<"\",\"boundary_times_p_numerator\":\""<<dec(lhs)<<"\",\"candidate_holds\":"<<(cmp(lhs,rhs)<=0?"true":"false");if(cross){U d=direct_h(S);std::cout<<",\"direct_hS\":\""<<dec(d)<<"\",\"cross_hS_equal\":"<<(cmp(d,hs)==0?"true":"false");}std::cout<<"}\n";return 0;
 }catch(const std::exception&e){std::cerr<<"completion_ie_selected: "<<e.what()<<"\n";return 1;}}
