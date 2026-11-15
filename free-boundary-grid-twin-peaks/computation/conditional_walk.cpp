// Symmetric transposition proposals rejected unless exactly two peaks remain.
// Uniform conditional measure is stationary; mixing/irreducibility is NOT asserted.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
int main(int argc,char**argv){
 if(argc<6||argc>8)return 2;int n=std::stoi(argv[1]),N=n*n;uint64_t steps=std::stoull(argv[2]),seed=std::stoull(argv[3]);std::string kind=argv[4];int stride=std::stoi(argv[5]);
 int reflect_every=argc>=7?std::stoi(argv[6]):0;
 int top_block_every=argc==8?std::stoi(argv[7]):0;
 if(n<3||steps<100||stride<1||reflect_every<0||top_block_every<0)return 2;
 std::mt19937_64 rng(seed);std::vector<std::vector<int>>adj(N);std::vector<std::pair<int,int>>edges;
 for(int v=0;v<N;v++)for(int u: {v-1,v+1,v-n,v+n})if(u>=0&&u<N&&abs(v%n-u%n)+abs(v/n-u/n)==1){adj[v].push_back(u);if(v<u)edges.emplace_back(v,u);}
 std::vector<int>order{0,kind=="far"?N-1:2},used(N),labels(N),loc(N+1);used[order[0]]=used[order[1]]=1;
 while((int)order.size()<N){std::vector<int>front;for(int v=0;v<N;v++)if(!used[v]&&std::any_of(adj[v].begin(),adj[v].end(),[&](int u){return used[u];}))front.push_back(v);int v=front[rng()%front.size()];used[v]=1;order.push_back(v);}
 for(int i=0;i<N;i++){labels[order[i]]=N-i;loc[N-i]=order[i];}
 auto peak=[&](int v){for(int u:adj[v])if(labels[u]>labels[v])return false;return true;};
 std::vector<int>pk(N);for(int v=0;v<N;v++)pk[v]=peak(v);
 uint64_t accepted=0,samples=0,reflection_attempts=0,reflection_accepts=0,block_attempts=0,block_accepts=0,top_moves=0;std::vector<uint64_t>hist(2*n),half1(2*n),half2(2*n);long double sqsum=0;int maxdist=0;std::vector<int>best=labels;
 for(uint64_t it=0;it<steps;it++){
  bool blockmove=top_block_every>0&&rng()%top_block_every==0;
  bool reflected=!blockmove&&reflect_every>0&&rng()%reflect_every==0;
  if(blockmove){
   block_attempts++;int oldtop=loc[N],x=oldtop%n,y=oldtop/n;
   std::vector<int>xs,ys;if(x>0)xs.push_back(x-1);if(x<n-1)xs.push_back(x);if(y>0)ys.push_back(y-1);if(y<n-1)ys.push_back(y);
   int bx=xs[rng()%xs.size()],by=ys[rng()%ys.size()];
   std::array<int,4>vs{by*n+bx,by*n+bx+1,(by+1)*n+bx,(by+1)*n+bx+1},oldval{},newval{};
   for(int i=0;i<4;i++)oldval[i]=labels[vs[i]];newval=oldval;std::shuffle(newval.begin(),newval.end(),rng);
   std::vector<int>affected;auto add=[&](int v){if(std::find(affected.begin(),affected.end(),v)==affected.end())affected.push_back(v);};
   for(int v:vs){add(v);for(int u:adj[v])add(u);}
   int old=0;for(int v:affected)old+=pk[v];int newtop=oldtop;
   for(int i=0;i<4;i++){labels[vs[i]]=newval[i];if(newval[i]==N)newtop=vs[i];}
   std::vector<int>npk;int now=0;for(int v:affected){npk.push_back(peak(v));now+=npk.back();}
   auto blocks=[&](int v){return (1+(v%n>0&&v%n<n-1))*(1+(v/n>0&&v/n<n-1));};
   double ratio=double(blocks(oldtop))/blocks(newtop);
   if(old==now&&(ratio>=1||std::generate_canonical<double,53>(rng)<ratio)){
    accepted++;block_accepts++;top_moves+=oldtop!=newtop;
    for(size_t i=0;i<affected.size();i++)pk[affected[i]]=npk[i];for(int v:vs)loc[labels[v]]=v;
   }else for(int i=0;i<4;i++)labels[vs[i]]=oldval[i];
  }else
  if(reflected){
   reflection_attempts++;int axis=rng()%2,l=rng()%n,r=rng()%n;if(l>r)std::swap(l,r);
   auto apply=[&](){for(int a=l;a<=(l+r)/2;a++)for(int b=0;b<n;b++){
    int v=axis==0?b*n+a:a*n+b,u=axis==0?b*n+(l+r-a):(l+r-a)*n+b;std::swap(labels[u],labels[v]);}};
   apply();std::vector<int>npk(N);int count=0;for(int w=0;w<N;w++)count+=(npk[w]=peak(w));
   if(count==2){pk=npk;reflection_accepts++;accepted++;for(int w=0;w<N;w++)loc[labels[w]]=w;}
   else apply();
  }else{
  int u,v;int mode=rng()%4;
  if(mode==0){int r=1+rng()%(N-1);u=loc[r];v=loc[r+1];}
  else if(mode==1){auto e=edges[rng()%edges.size()];u=e.first;v=e.second;}
  else {u=rng()%N;v=rng()%N;}
  std::array<int,10>affected{};int na=0;
  auto add=[&](int x){for(int i=0;i<na;i++)if(affected[i]==x)return;affected[na++]=x;};
  add(u);add(v);for(int w:adj[u])add(w);for(int w:adj[v])add(w);
  int old=0;for(int i=0;i<na;i++)old+=pk[affected[i]];
  std::swap(labels[u],labels[v]);std::array<int,10>newpk{};int now=0;
  for(int i=0;i<na;i++)now+=(newpk[i]=peak(affected[i]));
  if(old==now){accepted++;for(int i=0;i<na;i++)pk[affected[i]]=newpk[i];loc[labels[u]]=u;loc[labels[v]]=v;}
  else std::swap(labels[u],labels[v]);
  }
  if((it+1)%100000==0){int K=0;for(int w=0;w<N;w++){bool truth=peak(w);if(truth!=bool(pk[w]))return 4;K+=truth;}if(K!=2)return 5;}
  if(it>=steps/5&&(it-steps/5)%stride==0){int p=-1,q=-1;for(int w=0;w<N;w++)if(pk[w]){if(p<0)p=w;else q=w;}
   int dx=abs(p%n-q%n),dy=abs(p/n-q/n),d=dx+dy;hist[d]++;(it<steps*3/5?half1:half2)[d]++;sqsum+=dx*dx+dy*dy;samples++;
   if(d>maxdist){maxdist=d;best=labels;}
  }
 }
 auto array=[&](const auto&a){std::cout<<"[";for(size_t i=0;i<a.size();i++){if(i)std::cout<<",";std::cout<<a[i];}std::cout<<"]";};
 std::cout<<"{\"n\":"<<n<<",\"steps\":"<<steps<<",\"seed\":"<<seed<<",\"start\":\""<<kind<<"\",\"stride\":"<<stride<<",\"top_block_every\":"<<top_block_every<<",\"block_attempts\":"<<block_attempts<<",\"block_accepts\":"<<block_accepts<<",\"top_moves\":"<<top_moves<<",\"reflection_every\":"<<reflect_every<<",\"reflection_attempts\":"<<reflection_attempts<<",\"reflection_accepts\":"<<reflection_accepts<<",\"accepted\":"<<accepted<<",\"samples\":"<<samples<<",\"Q_empirical_unverified_mixing\":"<<double(sqsum/(N*samples))<<",\"distance_histogram\":";array(hist);std::cout<<",\"first_half_histogram\":";array(half1);std::cout<<",\"second_half_histogram\":";array(half2);std::cout<<",\"best_distance\":"<<maxdist<<",\"best_labels\":";array(best);std::cout<<"}\n";
}
