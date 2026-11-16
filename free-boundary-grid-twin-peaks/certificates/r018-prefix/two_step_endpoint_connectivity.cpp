// Exact communication census for endpoint-checked Q^2 on valid prefix words.
// Intermediate words are arbitrary injective words; only the two endpoints are
// required to satisfy the prescribed birth ranks.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using Word = std::vector<int>;
struct DSU { std::vector<int> p,s; explicit DSU(int n):p(n),s(n,1){std::iota(p.begin(),p.end(),0);} int f(int x){while(x!=p[x]){p[x]=p[p[x]];x=p[x];}return x;} void u(int a,int b){a=f(a);b=f(b);if(a==b)return;if(s[a]<s[b])std::swap(a,b);p[b]=a;s[a]+=s[b];} };
struct Grid { int n,N; std::vector<uint32_t>a; explicit Grid(int nn):n(nn),N(nn*nn),a(N){for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto [dx,dy]:{std::pair{1,0},{-1,0},{0,1},{0,-1}}){int X=x+dx,Y=y+dy;if(0<=X&&X<n&&0<=Y&&Y<n)a[v]|=uint32_t(1)<<(Y*n+X);}}} };
static uint64_t enc(const Word&w){uint64_t z=0;for(size_t i=0;i<w.size();i++)z|=uint64_t(w[i])<<(4*i);return z;}
static Word dec(uint64_t z,int r){Word w(r);for(int i=0;i<r;i++)w[i]=(z>>(4*i))&15;return w;}
static bool valid(const Word&w,int tau,const Grid&g){uint32_t seen=0;for(int t=1;t<=int(w.size());t++){int v=w[t-1];if(seen&(uint32_t(1)<<v))return false;bool b=(g.a[v]&seen)==0;if(b!=(t==1||t==tau))return false;seen|=uint32_t(1)<<v;}return true;}
static Word proposal(Word w,int i,int v){auto it=std::find(w.begin(),w.end(),v);if(it==w.end())w[i]=v;else std::swap(w[i],*it);return w;}
static std::string wordstr(const Word&w,const Grid&g){std::string s="[";for(size_t i=0;i<w.size();i++){if(i)s+=",";s+="["+std::to_string(w[i]%g.n)+","+std::to_string(w[i]/g.n)+"]";}return s+"]";}
struct Row{int n,r,tau,states,components;std::map<int,int>spectrum;std::vector<std::string>reps;};
static Row one(const Grid&g,int r,int tau){
 std::vector<uint64_t> st;Word w;
 auto rec=[&](auto&&self,uint32_t seen)->void{int t=int(w.size())+1;if(t>r){st.push_back(enc(w));return;}bool req=t==1||t==tau;for(int v=0;v<g.N;v++)if(!(seen&(uint32_t(1)<<v))&&(((g.a[v]&seen)==0)==req)){w.push_back(v);self(self,seen|(uint32_t(1)<<v));w.pop_back();}};rec(rec,0);
 std::sort(st.begin(),st.end());DSU d(int(st.size()));
 // For every Q-neighbour y of a valid endpoint x, record one endpoint owner.
 // Reencountering y realizes x --Q--> y --Q--> z and joins x,z.
 std::unordered_map<uint64_t,int> owner;owner.reserve(st.size()*size_t(r)*g.N*3/2+1);
 for(int id=0;id<int(st.size());id++){Word x=dec(st[id],r);for(int i=0;i<r;i++)for(int v=0;v<g.N;v++){uint64_t y=enc(proposal(x,i,v));auto [it,newkey]=owner.emplace(y,id);if(!newkey)d.u(id,it->second);}}
 std::map<int,std::vector<int>> groups;for(int i=0;i<int(st.size());i++)groups[d.f(i)].push_back(i);
 Row out{g.n,r,tau,int(st.size()),int(groups.size())};for(auto&[q,m]:groups){out.spectrum[int(m.size())]++;if(out.components>1)out.reps.push_back(wordstr(dec(st[m[0]],r),g));}return out;
}
static void quote(std::ostream&o,const std::string&s){o<<'"';for(char c:s){if(c=='"'||c=='\\')o<<'\\';o<<c;}o<<'"';}
static void emit(std::ostream&o,const Row&z){o<<"{\"n\":"<<z.n<<",\"r\":"<<z.r<<",\"tau\":"<<z.tau<<",\"states\":"<<z.states<<",\"components\":"<<z.components<<",\"size_spectrum\":{";bool f=true;for(auto[a,b]:z.spectrum){if(!f)o<<",";f=false;quote(o,std::to_string(a));o<<":"<<b;}o<<"},\"representatives\":[";for(size_t i=0;i<z.reps.size();i++){if(i)o<<",";quote(o,z.reps[i]);}o<<"]}";}
int main(){std::vector<Row> rows;for(int n:{2,3}){Grid g(n);for(int r=1;r<=g.N;r++)for(int tau=2;tau<=g.N;tau++)rows.push_back(one(g,r,tau));}std::cout<<"{\"method\":\"exact valid endpoints; union endpoints sharing an arbitrary injective one-step Q intermediate\",\"rows\":[";for(size_t i=0;i<rows.size();i++){if(i)std::cout<<",";emit(std::cout,rows[i]);}std::cout<<"],\"first_nontrivial_disconnected\":";bool found=false;for(const auto&z:rows)if(z.states>1&&z.components>1){emit(std::cout,z);found=true;break;}if(!found)std::cout<<"null";std::cout<<"}\n";}
