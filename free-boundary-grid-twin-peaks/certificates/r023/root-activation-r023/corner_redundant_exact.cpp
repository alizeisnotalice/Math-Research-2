// Restriction of the frozen all-state scanner to y=(0,0), x=(1,1) in C.
// Reuses its exact U128 recurrence and U512 fraction comparisons unchanged.
#define main frozen_all_state_main
#include "root_activation_exact.cpp"
#undef main

int main(int argc,char**argv) {
  if(argc!=3) return 2; int n=std::atoi(argv[1]); if(n<3||n>5)return 2;
  Grid g(n); const uint64_t cap=uint64_t(g.full)+1;
  std::vector<U128> h(cap); h[g.full]=1;
  for(uint32_t s=g.full-1;s;s--) { uint32_t b=g.nbr(s)&~s; U128 z=0;
    while(b){uint32_t bit=b&-b;b^=bit;z+=h[s|bit];}h[s]=z; }
  const int y=0,x=n+1; const uint32_t xbit=uint32_t(1)<<x;
  const uint32_t avail=g.full^((uint32_t(1)<<y)|g.adj[y]);
  std::vector<Rec> mp(g.N+1),ml(g.N+1),ms(g.N+1);std::vector<uint64_t> cnt(g.N+1);uint64_t checked=0;
  for(uint32_t C=avail;C;C=(C-1)&avail) {
    if(!(C&xbit)||!g.connected(C))continue;
    int k=__builtin_popcount(C);uint32_t S=C|(uint32_t(1)<<y);U128 num=0;uint32_t zset=g.adj[y];
    while(zset){uint32_t bit=zset&-zset;zset^=bit;num+=h[S|bit];}
    Rec r{true,C,y,num,h[S]};update(mp[k],r,k,0);update(ml[k],r,k,1);update(ms[k],r,k,2);cnt[k]++;checked++;
  }
  std::ofstream o(argv[2]);o<<"{\"n\":"<<n<<",\"fixed_y\":[0,0],\"forced_x\":[1,1]"
    <<",\"method\":\"exact corner-redundant subclass of the frozen U128 scanner\""
    <<",\"eligible_pairs\":"<<checked<<",\"by_C_size\":[";bool first=true;
  for(int k=1;k<=g.N;k++)if(cnt[k]){if(!first)o<<',';first=false;o<<"{\"C_size\":"<<k<<",\"eligible_pairs\":"<<cnt[k]
    <<",\"minimum_p_act\":";emit_rec(o,mp[k],k,g,h);o<<",\"minimum_p_times_C_plus_1\":";emit_rec(o,ml[k],k,g,h);
    o<<",\"minimum_p_times_sqrt_C\":";emit_rec(o,ms[k],k,g,h);o<<'}';}
  o<<"]}\n";
}
