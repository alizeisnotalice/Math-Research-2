#include <cassert>
#include <cstdint>
#include <vector>
#include <algorithm>
using U=uint64_t;
struct S{std::vector<U>x;size_t i=0;U operator()(){return x.at(i++);}};
static U below(S&s,U m){U cut=(U(0)-m)%m;for(;;){U x=s();if(x>=cut)return x%m;}}
int main(){
 {S s{{0}};assert(below(s,1)==0&&s.i==1);}
 {S s{{0,1}};assert(below(s,3)==1&&s.i==2);} // cut=1: raw 0 rejected
 {U m=(U(1)<<63)+1,cut=(U(0)-m)%m;assert(cut==(U(1)<<63)-1);S s{{cut-1,cut,cut+1}};assert(below(s,m)==cut&&s.i==2);assert(below(s,m)==cut+1&&s.i==3);}
 {std::vector<U> c{0,2,2,5};for(U u:{0,1,2,3,4}){int j=int(std::upper_bound(c.begin(),c.end(),u)-c.begin());assert(j==(u<2?1:3));}}
}
