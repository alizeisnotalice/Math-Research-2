#include <metal_stdlib>
using namespace metal;

// Development-only n<=16 move kernel.  The host supplies p.r as the current
// completed prefix length and never dispatches this at the terminal rank.
constant uint WORDS=8, MAXN=256;
struct P { uint n,N,tau,r,lo,hi; };
struct S { uint occ[WORDS], bd[WORDS]; ushort x[MAXN]; uint root1,root2; };
inline ulong seed(constant P&p){return(ulong(p.hi)<<32)|p.lo;}
inline ulong next64(thread ulong&x){x+=0x9e3779b97f4a7c15ul;ulong z=x;z=(z^(z>>30))*0xbf58476d1ce4e5b9ul;z=(z^(z>>27))*0x94d049bb133111ebul;return z^(z>>31);}
inline ulong below(thread ulong&x,ulong m){ulong t=(0ul-m)%m;while(1){ulong z=next64(x);if(z>=t)return z%m;}}
inline ulong moveKey(constant P&p,uint child,uint step){return seed(p)^(ulong(p.r)<<32)^(ulong(child)*0x9e3779b97f4a7c15ul)^(ulong(step+3u)*0xd1b54a32d192ed03ul);}
inline bool has(const thread S&s,uint v){return (s.occ[v>>5]&(1u<<(v&31)))!=0;}
inline bool inbd(const thread S&s,uint v){return (s.bd[v>>5]&(1u<<(v&31)))!=0;}
inline void put(thread uint*a,uint v){a[v>>5]|=1u<<(v&31);} inline void del(thread uint*a,uint v){a[v>>5]&=~(1u<<(v&31));}
inline void add(thread S&s,uint v,constant P&p){put(s.occ,v);del(s.bd,v);uint X=v%p.n,Y=v/p.n;if(X&&!has(s,v-1))put(s.bd,v-1);if(X+1<p.n&&!has(s,v+1))put(s.bd,v+1);if(Y&&!has(s,v-p.n))put(s.bd,v-p.n);if(Y+1<p.n&&!has(s,v+p.n))put(s.bd,v+p.n);}

// Canonical full ordered-word legality/rebuild.  It deliberately makes no
// unordered-subset shortcut: a swap is legal only after this position scan.
inline bool rebuild(thread S&s,constant P&p){
  for(uint w=0;w<WORDS;w++){s.occ[w]=0;s.bd[w]=0;}s.root1=0;s.root2=0;uint births=0;
  for(uint q=0;q<p.r;q++){
    uint v=uint(s.x[q]); if(v>=p.N||has(s,v))return false;
    bool birth=!inbd(s,v); if(birth){births++;if(q==0)s.root1=v;else if(q+1==p.tau)s.root2=v;}
    // Position-wise target; this also rejects an extra early/late birth.
    bool want=(q==0)||(q+1==p.tau);if(birth!=want)return false;add(s,v,p);
  }
  return births==((p.r<p.tau)?1u:2u);
}

// The SMC propagator supplies the already audited integer-CDF parent and its
// domain-2 allowed child vertex.  This performs real genealogy, rather than
// reconstructing an arbitrary ordering from an occupied subset.
kernel void orderedPrefixCopyAppend(device S*src[[buffer(0)]],device uint*parent[[buffer(1)]],device uint*childVertex[[buffer(2)]],device S*dst[[buffer(3)]],constant P&p[[buffer(4)]],constant uint&particles[[buffer(5)]],device atomic_uint*bad[[buffer(6)]],uint child[[thread_position_in_grid]]){
  if(child>=particles)return;thread S q=src[parent[child]];q.x[p.r-1]=ushort(childVertex[child]);
  // An audited propagator only supplies valid edges.  Keep this guard in the
  // development kernel so a prefix/order bug cannot silently enter a move.
  if(!rebuild(q,p)){atomic_store_explicit(bad,1u,memory_order_relaxed);dst[child]=src[parent[child]];return;}dst[child]=q;
}

kernel void orderedPrefixMove(device S*state[[buffer(0)]],constant P&p[[buffer(1)]],constant uint&particles[[buffer(2)]],constant uint&L[[buffer(3)]],uint child[[thread_position_in_grid]]){
  if(child>=particles||!L)return;thread S s=state[child];
  for(uint step=0;step<L;step++){
    ulong rng=moveKey(p,child,step);uint i=uint(below(rng,ulong(p.r))),v=uint(below(rng,ulong(p.N))),j=p.r;
    for(uint k=0;k<p.r;k++)if(uint(s.x[k])==v){j=k;break;}
    thread S trial=s;if(j<p.r){ushort t=trial.x[i];trial.x[i]=trial.x[j];trial.x[j]=t;}else trial.x[i]=ushort(v);
    if(rebuild(trial,p))s=trial;
  }
  state[child]=s;
}
