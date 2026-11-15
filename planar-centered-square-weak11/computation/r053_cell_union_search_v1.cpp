#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <cstring>
struct R{double x0,x1,y0,y1;};
static double ar(const double *w,const int *p){
 R r[9];double xs[18],s=0;int n=0;
 const double a[]={0,.25,.75},b[]={.25,.75,1};
 for(int k=0;k<9;k++){int j=p[k];if(w[j]<=0)continue;s+=w[j];double h=sqrt(s)/2;int x=j/3,y=j%3;r[n]={a[x]-h,b[x]+h,a[y]-h,b[y]+h};xs[2*n]=r[n].x0;xs[2*n+1]=r[n].x1;n++;}
 std::sort(xs,xs+2*n);double ans=0;
 for(int k=1;k<2*n;k++){if(xs[k]<=xs[k-1])continue;std::pair<double,double>v[9];int m=0;for(int i=0;i<n;i++)if(r[i].x0<=xs[k-1]&&r[i].x1>=xs[k])v[m++]={r[i].y0,r[i].y1};if(!m)continue;std::sort(v,v+m);double lo=v[0].first,hi=v[0].second,len=0;for(int i=1;i<m;i++){if(v[i].first>hi){len+=hi-lo;lo=v[i].first;hi=v[i].second;}else hi=std::max(hi,v[i].second);}len+=hi-lo;ans+=(xs[k]-xs[k-1])*len;}
 return ans;
}
extern "C" double area(const double*w,const int*p){return ar(w,p);}
extern "C" double full(const double*w,int*out,int maximize){int p[9];std::iota(p,p+9,0);double best=maximize?-1:100;do{double z=ar(w,p);if((maximize&&z>best)||(!maximize&&z<best)){best=z;std::copy(p,p+9,out);}}while(std::next_permutation(p,p+9));return best;}
extern "C" double heuristic(const double*w,int*out,int seeds,unsigned long seed){std::mt19937_64 rng(seed);int p[9];double best=100;
 for(int rep=0;rep<seeds;rep++){std::iota(p,p+9,0);if(rep==0)std::stable_sort(p,p+9,[&](int i,int j){return w[i]<w[j];});else std::shuffle(p,p+9,rng);double z=ar(w,p);bool improve=true;while(improve){improve=false;double v=z;int bi=-1,bj=-1;for(int i=0;i<9;i++)for(int j=i+1;j<9;j++){std::swap(p[i],p[j]);double q=ar(w,p);std::swap(p[i],p[j]);if(q<v-1e-13){v=q;bi=i;bj=j;}}if(bi>=0){std::swap(p[bi],p[bj]);z=v;improve=true;}}if(z<best){best=z;std::copy(p,p+9,out);}}
 return best;}
