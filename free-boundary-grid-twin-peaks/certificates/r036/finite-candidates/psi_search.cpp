// Low-memory one-peak Ψ search.  Search arithmetic is double; verifier is separate.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
struct State{int n,N,root;std::vector<int>a,pos,l,b,sole;double psi=0;std::vector<std::vector<int>> nb;
 State(int q,int r):n(q),N(q*q),root(r),a(N),pos(N),l(N),b(N),sole(N),nb(N){for(int v=0;v<N;v++){int x=v%n,y=v/n;for(auto d:{std::pair{-1,0},{1,0},{0,-1},{0,1}}){int X=x+d.first,Y=y+d.second;if(X>=0&&X<n&&Y>=0&&Y<n)nb[v].push_back(Y*n+X);}}}
 bool eligible(int p)const{return p!=root&&std::find(nb[root].begin(),nb[root].end(),p)==nb[root].end();}
 int earlycount(int v)const{int c=0;for(int u:nb[v])c+=pos[u]<pos[v];return c;}
 int soleparent(int v)const{int z=-1;for(int u:nb[v])if(pos[u]<pos[v]){if(z>=0)return -1;z=u;}return z;}
 void recompute_p(int p){l[p]=N+1;for(int u:nb[p])l[p]=std::min(l[p],pos[u]);b[p]=N+1;for(int v:nb[p])if(sole[v]==p)b[p]=std::min(b[p],pos[v]);}
 double term(int p)const{if(!eligible(p))return 0;int R=b[p]-l[p]-1;if(R<=0)throw std::runtime_error("bad R");return double(l[p]-1)/R;}
 void init_bfs(){std::vector<int>d(N,-1);d[root]=0;a[0]=root;for(int i=0,k=1;i<k;i++)for(int u:nb[a[i]])if(d[u]<0){d[u]=d[a[i]]+1;a[k++]=u;}for(int i=0;i<N;i++)pos[a[i]]=i+1;rebuild();}
 void init_row(){if(root!=0)throw std::runtime_error("row root");for(int i=0;i<N;i++)a[i]=i,pos[i]=i+1;rebuild();}
 void init_ring(){if(root!=0)throw std::runtime_error("ring root");int z=0;a[z++]=0;for(int k=1;k<n;k++){for(int y=0;y<=k;y++)a[z++]=y*n+k;for(int x=k-1;x>=0;x--)a[z++]=k*n+x;}for(int i=0;i<N;i++)pos[a[i]]=i+1;rebuild();}
 void rebuild(){for(int v=0;v<N;v++){if(v!=root&&earlycount(v)==0)throw std::runtime_error("not onepeak");sole[v]=soleparent(v);}psi=0;for(int p=0;p<N;p++){recompute_p(p);psi+=term(p);}}
 bool swap_try(int i,int j){if(i==j||i==0||j==0)return false;std::vector<int>S={a[i],a[j]};size_t ns=S.size();for(size_t z=0;z<ns;z++)for(int u:nb[S[z]])S.push_back(u);std::sort(S.begin(),S.end());S.erase(std::unique(S.begin(),S.end()),S.end());std::vector<int>P=S;for(int v:S)if(sole[v]>=0)P.push_back(sole[v]);size_t np=P.size();for(size_t z=0;z<np;z++)for(int u:nb[P[z]])P.push_back(u);std::sort(P.begin(),P.end());P.erase(std::unique(P.begin(),P.end()),P.end());double old=0;for(int p:P)old+=term(p);std::swap(a[i],a[j]);pos[a[i]]=i+1;pos[a[j]]=j+1;for(int v:S)sole[v]=soleparent(v);for(int v:S)if(v!=root&&earlycount(v)==0)return false;for(int p:P)recompute_p(p);double now=0;for(int p:P)now+=term(p);psi+=now-old;return true;}
};
int main(int c,char**v){try{if(c<6||c>8)throw std::runtime_error("usage n root seed proposals out [audit] [bfs|row|ring]");int n=std::stoi(v[1]),r=std::stoi(v[2]);uint64_t seed=std::stoull(v[3]),steps=std::stoull(v[4]);if(n<2||r<0||r>=n*n||steps==0)throw std::runtime_error("bounds");int chk=c>=7?std::stoi(v[6]):1000;if(chk<1)throw std::runtime_error("audit");State s(n,r);std::string init=c==8?v[7]:"bfs";if(init=="bfs")s.init_bfs();else if(init=="row")s.init_row();else if(init=="ring")s.init_ring();else throw std::runtime_error("init");auto audit=[](State& x){double old=x.psi;State z=x;z.rebuild();double tol=1e-11*std::max(1.0,std::abs(z.psi));if(std::abs(old-z.psi)>tol)throw std::runtime_error("incremental mismatch");x=std::move(z);};State best=s;std::mt19937_64 g(seed);uint64_t acc=0;for(uint64_t z=0;z<steps;z++){int i=1+g()%(s.N-1),j=1+g()%(s.N-1);double before=s.psi;State old=s;if(s.swap_try(i,j)){if(z%chk==0)audit(s);double after=s.psi;if(after<=before||std::log(std::generate_canonical<double,53>(g))<-(after-before)/.05)acc++;else s=std::move(old);}else s=std::move(old);if(s.psi<best.psi)best=s;}audit(s);audit(best);std::ofstream f(v[5]);if(!f)throw std::runtime_error("out");f.precision(17);f<<"n "<<n<<" root "<<r<<" seed "<<seed<<" proposals "<<steps<<" accepted "<<acc<<" psi "<<best.psi<<" norm "<<best.psi/(best.N*std::log(double(n)))<<" init "<<init<<"\nword";for(int x:best.a)f<<' '<<x;f<<'\n';}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
