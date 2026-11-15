#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>

int main(){
 constexpr int n=5;
 // Exact n=5 zero-orbit core from orbit-search-n5-d3.json.
 const int core[25]={25,7,24,4,16,9,13,23,20,22,21,6,12,15,19,3,1,14,17,5,10,8,11,2,18};
 std::array<int,n> px,py;std::iota(px.begin(),px.end(),0);
 uint64_t configs=0,k2=0;int min_peaks=99,max_k2_distance=0;bool fail=false;std::array<int,4> witness{};
 do{
  std::iota(py.begin(),py.end(),0);
  do{
   int a[25];for(int y=0;y<n;y++)for(int x=0;x<n;x++)a[y*n+x]=core[py[y]*n+px[x]];
   for(int mask=0;mask<256;mask++){
    int pp[25],count=0;
    for(int y=0;y<n;y++)for(int x=0;x<n;x++){
     bool peak=true;auto lower=[&](int X,int Y,int gap){if(X<0||X>=n||Y<0||Y>=n)return; if(mask&(1<<gap))return; if(a[y*n+x]<a[Y*n+X])peak=false;};
     if(x)lower(x-1,y,x-1);if(x+1<n)lower(x+1,y,x);if(y)lower(x,y-1,4+y-1);if(y+1<n)lower(x,y+1,4+y);
     if(peak)pp[count++]=y*n+x;
    }
    configs++;min_peaks=std::min(min_peaks,count);
    if(count==2){
      k2++;int x0=pp[0]%n,x1=pp[1]%n,y0=pp[0]/n,y1=pp[1]/n;bool between=false;
      max_k2_distance=std::max(max_k2_distance,std::abs(x0-x1)+std::abs(y0-y1));
      for(int g=std::min(x0,x1);g<std::max(x0,x1);g++)if(mask&(1<<g))between=true;
      for(int g=std::min(y0,y1);g<std::max(y0,y1);g++)if(mask&(1<<(4+g)))between=true;
      if(between&&!fail){fail=true;witness={mask,pp[0],pp[1],int(px[0]+5*py[0])};}
    }
   }
  }while(std::next_permutation(py.begin(),py.end()));
 }while(std::next_permutation(px.begin(),px.end()));
 std::cout<<"configs="<<configs<<" min_peaks="<<min_peaks<<" k2="<<k2<<" max_k2_distance="<<max_k2_distance<<" gap_between_k2="<<(fail?"true":"false");
 if(fail)std::cout<<" witness_mask="<<witness[0]<<" peaks="<<witness[1]<<","<<witness[2];
 std::cout<<"\n";
}
