#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>
struct Rect{double a,b,c,d;};
static double area(const std::vector<Rect>&r){std::vector<double>xs;for(auto q:r){if(q.a<q.b&&q.c<q.d){xs.push_back(q.a);xs.push_back(q.b);}}std::sort(xs.begin(),xs.end());double ans=0;for(size_t k=1;k<xs.size();k++){if(xs[k]<=xs[k-1])continue;std::vector<std::pair<double,double>>ys;for(auto q:r)if(q.a<=xs[k-1]&&q.b>=xs[k]&&q.c<q.d)ys.push_back({q.c,q.d});if(ys.empty())continue;std::sort(ys.begin(),ys.end());double l=ys[0].first,h=ys[0].second,t=0;for(auto y:ys){if(y.first>h){t+=h-l;l=y.first;h=y.second;}else h=std::max(h,y.second);}ans+=(xs[k]-xs[k-1])*(t+h-l);}return ans;}
int main(int argc,char**argv){if(argc!=4)return 2;std::ifstream in(argv[1]);std::ofstream out(argv[2]);int maxdel=std::stoi(argv[3]),id,n;long long Q;const double B=3.299999995647862;out<<std::setprecision(17);
 while(in>>id>>n>>Q){int count=1<<n,all=count-1;std::vector<long long>x(n),y(n),w(n);long long W=0;for(int j=0;j<n;j++){in>>x[j]>>y[j]>>w[j];W+=w[j];}std::vector<long long>mass(count),xmin(count),xmax(count),ymin(count),ymax(count);std::vector<int>eligible,rank(count),pop(count);std::vector<std::vector<long long>>cap(n,std::vector<long long>(count));
 for(int s=1;s<count;s++){int j=__builtin_ctz((unsigned)s),t=s^(1<<j);mass[s]=mass[t]+w[j];pop[s]=pop[t]+1;xmin[s]=t?std::min(xmin[t],x[j]):x[j];xmax[s]=t?std::max(xmax[t],x[j]):x[j];ymin[s]=t?std::min(ymin[t],y[j]):y[j];ymax[s]=t?std::max(ymax[t],y[j]):y[j];long long d=std::max(xmax[s]-xmin[s],ymax[s]-ymin[s]);if((__int128)d*d*W<(__int128)mass[s]*Q*Q){eligible.push_back(s);rank[s]=pop[s];for(int j=0;j<n;j++)if(s>>j&1)cap[j][s]=mass[s];}}
 for(int bit=0;bit<n;bit++)for(int s=0;s<count;s++)if(s>>bit&1){int t=s^(1<<bit);rank[s]=std::max(rank[s],rank[t]);for(int j=0;j<n;j++)cap[j][s]=std::max(cap[j][s],cap[j][t]);}
 auto square=[&](int j,long long m){double h=sqrt((double)m/W)/2;return Rect{(double)x[j]/Q-h,(double)x[j]/Q+h,(double)y[j]/Q-h,(double)y[j]/Q+h};};
 std::vector<double>leafU(count,-1),leafGarea(count);std::vector<std::vector<Rect>>leafG(count);std::vector<int>leafBound(count,-1);
 auto leaf=[&](int s){if(leafU[s]>=0)return;long long d=std::max(xmax[s]-xmin[s],ymax[s]-ymin[s]);double m=(double)mass[s]/W;double u=B*std::max(m,(double)d*d/Q/Q);int which=0;if(rank[s]*m<u){u=rank[s]*m;which=1;}std::vector<int>ids;for(int j=0;j<n;j++)if(s>>j&1){ids.push_back(j);leafG[s].push_back(square(j,w[j]));}if((__int128)d*d*W<(__int128)mass[s]*Q*Q){double h=sqrt(m)/2;leafG[s].push_back({(double)xmax[s]/Q-h,(double)xmin[s]/Q+h,(double)ymax[s]/Q-h,(double)ymin[s]/Q+h});}leafGarea[s]=area(leafG[s]);
 for(int type=0;type<4;type++){std::sort(ids.begin(),ids.end(),[&](int a,int b){if(type==0)return w[a]!=w[b]?w[a]<w[b]:a<b;if(type==1)return w[a]!=w[b]?w[a]>w[b]:a<b;if(type==2)return x[a]!=x[b]?x[a]<x[b]:a<b;return y[a]!=y[b]?y[a]<y[b]:a<b;});int prefix=0;std::vector<Rect>proxy;for(int j:ids){prefix|=1<<j;proxy.push_back(square(j,cap[j][prefix]));}double z=area(proxy);if(z<u){u=z;which=2+type;}}
 leafU[s]=u;leafBound[s]=which;};
 double best=1e100;int bestD=-1,checked=0;std::vector<int>bestLeaves;std::vector<double>bestParts;
 for(int D=0;D<count;D++){if(pop[D]>maxdel)continue;int rem=all^D,parent[16];std::iota(parent,parent+n,0);auto find=[&](int j){while(parent[j]!=j)j=parent[j];return j;};for(int e:eligible)if((e&rem)==e){int a=__builtin_ctz((unsigned)e);for(int j=a+1;j<n;j++)if(e>>j&1){int ra=find(a),rj=find(j);parent[ra]=rj;}}int masks[16]={};for(int j=0;j<n;j++)if(rem>>j&1)masks[find(j)]|=1<<j;std::vector<Rect>global;for(int j=0;j<n;j++)if(D>>j&1)global.push_back(square(j,cap[j][all]));double residual=0;std::vector<int>parts;std::vector<double>us;for(int s:masks)if(s){leaf(s);parts.push_back(s);us.push_back(leafU[s]);global.insert(global.end(),leafG[s].begin(),leafG[s].end());residual+=leafU[s]-leafGarea[s];}double z=area(global)+residual;checked++;if(z<best){best=z;bestD=D;bestLeaves=parts;bestParts=us;}}
 double actual=-1;if(best>3.5){std::vector<Rect>halo;for(int s:eligible){double h=sqrt((double)mass[s]/W)/2;halo.push_back({(double)xmax[s]/Q-h,(double)xmin[s]/Q+h,(double)ymax[s]/Q-h,(double)ymin[s]/Q+h});}actual=area(halo);}
 out<<id<<" "<<n<<" "<<best<<" "<<bestD<<" "<<rank[all]<<" "<<eligible.size()<<" "<<actual<<" "<<checked<<" "<<bestLeaves.size();for(size_t j=0;j<bestLeaves.size();j++)out<<" "<<bestLeaves[j]<<" "<<bestParts[j]<<" "<<leafBound[bestLeaves[j]];out<<"\n";
 }
}
