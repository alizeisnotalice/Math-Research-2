#define main production_main
#include "large_side_orbit_dp.cpp"
#undef main
// Independent lexicographic rank via binomial hockey-stick subtraction.
static U64 rank_independent(const std::vector<int>& x,int m){
 int l=static_cast<int>(x.size()), top=m+l-1, previous=-1; U64 rank=0;
 for(int i=0;i<l;++i){int y=x[i]+i, remaining=l-i;
  rank += binom(top-previous-1,remaining)-binom(top-y,remaining);
  previous=y;
 }
 return rank;
}
int main(){
 for(int l:{7,6}){
  int m=1<<(11-l); U64 total=binom(m+l-1,l),sum=0;std::vector<int> previous;
  for(U64 k=0;k<total;++k){auto x=unrank(m,l,k);
   if(!std::is_sorted(x.begin(),x.end()) || x.front()<0 || x.back()>=m || rank_independent(x,m)!=k || (!previous.empty()&&!(previous<x)))return 3;
   sum+=weight(x);previous=x;
  }
  if(sum!=(1ULL<<(l*(11-l))))return 4;
  std::cout<<"{\"sorted_side\":"<<l<<",\"representatives\":"<<total<<",\"weighted_graphs\":"<<sum<<",\"all_inverse_ranks_pass\":true}\n";
 }
}
