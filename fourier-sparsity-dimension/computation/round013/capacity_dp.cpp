// Exact support-partition capacity DP; states are codeword parameter masks.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std; using U=uint64_t;
int main(int ac,char**av){if(ac<3){cerr<<"usage: capacity_dp D basis_masks...\n";return 2;}int d=atoi(av[1]);if(d<1||d>20||ac!=d+2){cerr<<"bad input\n";return 2;}int M=1<<d;vector<U> cw(M);for(int i=0;i<d;i++){U g=strtoull(av[i+2],0,0);for(int s=0;s<(1<<i);s++)cw[s|(1<<i)]=cw[s]^g;}auto u=cw;sort(u.begin(),u.end());if(unique(u.begin(),u.end())!=u.end()){cerr<<"dependent basis\n";return 3;}vector<int> sub(M);for(int s=0;s<M;s++){int r=0;for(int a=0;a<M;a++)if((cw[a]&~cw[s])==0)r++;if(r<=0||(r&(r-1))){cerr<<"bad subcount\n";return 4;}sub[s]=__builtin_ctz((unsigned)r);}
 vector<int> order(M);iota(order.begin(),order.end(),0);sort(order.begin(),order.end(),[&](int x,int y){return __builtin_popcountll(cw[x])<__builtin_popcountll(cw[y]);});vector<int> dp(M,-1),pick(M);dp[0]=0;for(int s:order)if(s){U R=cw[s];int lo=__builtin_ctzll(R);for(int a=1;a<M;a++)if((cw[a]&~R)==0 && ((cw[a]>>lo)&1)){int v=sub[a]+1+dp[s^a];if(v>dp[s])dp[s]=v,pick[s]=a;}}
 int root=find(cw.begin(),cw.end(),~U(0))-cw.begin();if(root==M||dp[root]<0){cerr<<"full/root missing\n";return 5;}cout<<"{\"D\":"<<d<<",\"kappa\":"<<dp[root]<<",\"partition\":[";bool f=1;for(int s=root;s;){if(pick[s]<=0){cerr<<"invalid recursion\n";return 6;}if(!f)cout<<',';f=0;cout<<"\"0x"<<hex<<cw[pick[s]]<<dec<<"\"";s^=pick[s];}cout<<"]}\n";
}
