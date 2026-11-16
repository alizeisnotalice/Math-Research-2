import Foundation
import Metal

// Experimental original-law SMC, n<=32 (MAXN prefix capacity is 1024).  It stores an exact occupancy and
// external-boundary bitset.  It is deliberately separate from r013 production.
let metal = """
  #include <metal_stdlib>
  using namespace metal;
  constant uint W=32, MAXN=1024;
  struct P{uint n,N,tau,r,lo,hi;};
  struct S{uint o[W],b[W],r1,r2;};
  inline ulong seed(constant P&p){return(ulong(p.hi)<<32)|p.lo;}
  inline ulong nxt(thread ulong&x){x+=0x9e3779b97f4a7c15ul;ulong z=x;z=(z^(z>>30))*0xbf58476d1ce4e5b9ul;z=(z^(z>>27))*0x94d049bb133111ebul;return z^(z>>31);}
  inline ulong key(constant P&p,uint c,uint d){return seed(p)^(ulong(p.r)<<32)^(ulong(c)*0x9e3779b97f4a7c15ul)^(ulong(d)*0xd1b54a32d192ed03ul);}
  inline ulong below(thread ulong&x,ulong n){ulong t=(0ul-n)%n;while(1){ulong z=nxt(x);if(z>=t)return z%n;}}
  inline uint lm(constant P&p){uint r=p.N&31;return r?((1u<<r)-1u):0xffffffffu;}
  inline bool has(const thread S&s,uint v){return(s.o[v>>5]&(1u<<(v&31)))!=0;}
  inline void put(thread uint*x,uint v){x[v>>5]|=1u<<(v&31);} inline void del(thread uint*x,uint v){x[v>>5]&=~(1u<<(v&31));}
  inline uint aw(const thread S&s,uint q,constant P&p){uint last=(p.N-1)>>5;if(q>last)return 0;uint f=~s.o[q];if(q==last)f&=lm(p);if(p.r==1)return f;if(p.r==p.tau)return f&~s.b[q];return s.b[q];}
  inline uint cnt(const thread S&s,constant P&p){uint z=0;for(uint q=0;q<W;q++)z+=popcount(aw(s,q,p));return z;}
  inline uint kth(const thread S&s,uint k,constant P&p){for(uint q=0;q<W;q++){uint x=aw(s,q,p),c=popcount(x);if(k>=c){k-=c;continue;}for(uint j=0;j<32;j++)if(x&(1u<<j)){if(!k)return(q<<5)+j;--k;}}return p.N;}
  inline void add(thread S&s,uint v,constant P&p){put(s.o,v);del(s.b,v);uint x=v%p.n,y=v/p.n;if(x&&!has(s,v-1))put(s.b,v-1);if(x+1<p.n&&!has(s,v+1))put(s.b,v+1);if(y&&!has(s,v-p.n))put(s.b,v-p.n);if(y+1<p.n&&!has(s,v+p.n))put(s.b,v+p.n);}
  kernel void count(device S*x[[buffer(0)]],device uint*a[[buffer(1)]],constant P&p[[buffer(2)]],device uint*d[[buffer(3)]],uint i[[thread_position_in_grid]]){thread S q=x[i];a[i]=d[0]?0:cnt(q,p);}
  kernel void sum(device uint*a[[buffer(0)]],device ulong*g[[buffer(1)]],constant uint&n[[buffer(2)]],uint j[[thread_position_in_grid]]){uint lo=256*j,hi=min(lo+256,n);ulong z=0;for(uint i=lo;i<hi;i++)z+=a[i];g[j]=z;}
  kernel void cdf(device ulong*g[[buffer(0)]],device ulong*c[[buffer(1)]],device ulong*t[[buffer(2)]],device uint*d[[buffer(3)]],constant uint&G[[buffer(4)]],constant uint&r[[buffer(5)]],uint i[[thread_position_in_grid]]){if(i)return;ulong z=0;for(uint j=0;j<G;j++){z+=g[j];c[j]=z;}t[r]=z;if(!z)d[0]=1;}
  kernel void prop(device S*x[[buffer(0)]],device uint*a[[buffer(1)]],device ulong*c[[buffer(2)]],device ulong*t[[buffer(3)]],device S*y[[buffer(4)]],constant P&p[[buffer(5)]],constant uint&n[[buffer(6)]],constant uint&r[[buffer(7)]],device uint*d[[buffer(8)]],constant uint&st[[buffer(9)]],device uint*parents[[buffer(10)]],device ushort*px[[buffer(11)]],device ushort*py[[buffer(12)]],uint ch[[thread_position_in_grid]]){if(ch>=n)return;if(d[0]){parents[ch]=0xffffffffu;y[ch]=x[ch];for(uint k=0;k+1<p.r;k++)py[ch*MAXN+k]=px[ch*MAXN+k];return;}ulong T=t[r],z=key(p,ch,1),u=below(z,T);if(st)u=(ulong(ch)*T+u)/n;uint G=(n+255)/256,L=0,H=G-1;while(L<H){uint m=(L+H)/2;if(u<c[m])H=m;else L=m+1;}ulong pre=L?c[L-1]:0,v=u-pre,s=0;uint lo=256*L,hi=min(lo+256,n),par=lo;for(uint i=lo;i<hi;i++){if(v<s+a[i]){par=i;break;}s+=a[i];}parents[ch]=par;thread S q=x[par];for(uint k=0;k+1<p.r;k++)py[ch*MAXN+k]=px[par*MAXN+k];z=key(p,ch,2);uint w=kth(q,uint(below(z,a[par])),p);py[ch*MAXN+p.r-1]=ushort(w);if(p.r==1)q.r1=w;if(p.r==p.tau)q.r2=w;add(q,w,p);y[ch]=q;}
  inline bool rebuildPrefix(thread S&s,device ushort*pref,uint id,constant P&p){for(uint w=0;w<W;w++){s.o[w]=0;s.b[w]=0;}s.r1=0;s.r2=0;uint births=0;for(uint q=0;q<p.r;q++){uint v=uint(pref[id*MAXN+q]);if(v>=p.N||has(s,v))return false;bool birth=(s.b[v>>5]&(1u<<(v&31)))==0;bool want=(q==0)||(q+1==p.tau);if(birth!=want)return false;if(birth){births++;if(q==0)s.r1=v;else s.r2=v;}add(s,v,p);}return births==((p.r<p.tau)?1u:2u);}
  kernel void move(device S*state[[buffer(0)]],device ushort*pref[[buffer(1)]],constant P&p[[buffer(2)]],constant uint&n[[buffer(3)]],constant uint&L[[buffer(4)]],device uint*d[[buffer(5)]],uint ch[[thread_position_in_grid]]){if(ch>=n||!L||d[0])return;thread S s=state[ch];for(uint step=0;step<L;step++){ulong z=key(p,ch,step+3u);uint i=uint(below(z,ulong(p.r))),v=uint(below(z,ulong(p.N))),j=p.r;for(uint k=0;k<p.r;k++)if(uint(pref[ch*MAXN+k])==v){j=k;break;}ushort oldi=pref[ch*MAXN+i],oldj=(j<p.r?pref[ch*MAXN+j]:ushort(v));pref[ch*MAXN+i]=oldj;if(j<p.r)pref[ch*MAXN+j]=oldi;thread S trial=s;if(rebuildPrefix(trial,pref,ch,p))s=trial;else{pref[ch*MAXN+i]=oldi;if(j<p.r)pref[ch*MAXN+j]=oldj;}}state[ch]=s;}
  kernel void term(device S*x[[buffer(0)]],device uint*a[[buffer(1)]],device ulong*m[[buffer(2)]],device ulong*f[[buffer(3)]],device long*h[[buffer(4)]],constant P&p[[buffer(5)]],constant uint&n[[buffer(6)]],device uint*d[[buffer(7)]],uint j[[thread_position_in_grid]]){if(d[0]){m[j]=0;f[j]=0;h[j]=0;return;}uint lo=256*j,hi=min(lo+256,n);ulong z=0,ff=0;long hh=0;for(uint i=lo;i<hi;i++){thread S q=x[i];uint r=q.r2;if(p.r==p.tau)r=kth(q,0,p);int X=int(q.r1%p.n)-int(r%p.n),Y=int(q.r1/p.n)-int(r/p.n);z+=ulong(a[i])*ulong(X*X+Y*Y);if(p.tau<=p.N){int xp=2*int(q.r1%p.n)-int(p.n-1),yp=2*int(q.r1/p.n)-int(p.n-1),xq=2*int(r%p.n)-int(p.n-1),yq=2*int(r/p.n)-int(p.n-1);ff+=ulong(a[i])*ulong(xp*xp+yp*yp+xq*xq+yq*yq);hh+=long(a[i])*long(xp*xq+yp*yq);}}m[j]=z;f[j]=ff;h[j]=hh;}
  """
// A deliberately separate runtime library: terminal-only exact inverse fibres.
// Keeping this small avoids recompiling or perturbing the frozen SMC kernels.
let inverseMetal = """
  #include <metal_stdlib>
  using namespace metal;
  constant uint W2=32, MAXN2=1024;
  struct P2{uint n,N,tau,r,lo,hi;};
  struct S2{uint o[W2],b[W2],r1,r2;};
  inline uint md2(uint a,uint b,uint n){return uint(abs(int(a%n)-int(b%n))+abs(int(a/n)-int(b/n)));}
  inline bool ed2(uint a,uint b,uint n){return md2(a,b,n)==1;}
  inline bool has2(const thread S2&s,uint v){return (s.o[v>>5]&(1u<<(v&31)))!=0;}
  inline uint sv2(device const ushort*pref,const thread S2&q,uint id,uint k,uint N){
    if(k+1<N)return uint(pref[id*MAXN2+k]);
    for(uint v=0;v<N;v++)if(!has2(q,v))return v;
    return N;
  }
  // The rank table is terminal-only: rank[id*1024+v] is sigma^{-1}(v).
  // It converts every "earlier tau neighbour" query from an O(N) prefix scan
  // to inspection of at most the four graph neighbours of w.
  inline bool nearFast(device const ushort*rank,uint id,uint d,uint mask,uint r,uint w,uint N,uint n){
    uint ti=r-popcount(mask), x=w%n, y=w/n;
    uint nb[4]={x?w-1:N,x+1<n?w+1:N,y?w-n:N,y+1<n?w+n:N};
    for(uint j=0;j<4;j++){uint z=nb[j];if(z>=N)continue;uint rz=uint(rank[id*MAXN2+z]);
      if(rz==0 || (rz>=d && rz-d+1<ti) || (rz>0 && rz<d && (mask&(1u<<(rz-1)))))return true;
    }
    return false;
  }
  inline bool canon2(device const ushort*pref,const thread S2&q,uint id,uint N,uint n,uint d,uint v){
    uint cur=sv2(pref,q,id,0,N);
    for(uint k=1;k<d;k++){
      uint best=N,old=md2(cur,v,n);
      for(uint w=0;w<N;w++)if(ed2(cur,w,n)&&md2(w,v,n)+1==old&&w<best)best=w;
      if(best!=sv2(pref,q,id,k,N))return false;
      cur=best;
    }
    return true;
  }
  inline ulong invFast(device const ushort*pref,device const ushort*rank,const thread S2&q,uint id,uint N,uint n,uint d){
    if(md2(sv2(pref,q,id,0,N),sv2(pref,q,id,d-1,N),n)!=d-1)return 0;
    for(uint k=0;k+1<d;k++)if(!ed2(sv2(pref,q,id,k,N),sv2(pref,q,id,k+1,N),n))return 0;
    uint plen=d-1,full=(1u<<plen)-1,T=N-plen,u=sv2(pref,q,id,0,N); ulong ans=0;
    uint last=sv2(pref,q,id,d-1,N),x=last%n,y=last/n;
    uint cand[4]={x?last-1:N,x+1<n?last+1:N,y?last-n:N,y+1<n?last+n:N};
    for(uint ci=0;ci<4;ci++){uint v=cand[ci];
      if(v>=N||uint(rank[id*MAXN2+v])<d||md2(u,v,n)!=d||!canon2(pref,q,id,N,n,d,v))continue;
      ulong cur[4]={0,0,0,0},nxt[4];cur[0]=1;
      for(uint r=1;r<N;r++){
        for(uint z=0;z<4;z++)nxt[z]=0;
        for(uint mask=0;mask<=full;mask++){
          ulong c=cur[mask];if(!c)continue;uint ti=r-popcount(mask);
          for(uint j=0;j<plen;j++)if(!(mask&(1u<<j))&&nearFast(rank,id,d,mask,r,sv2(pref,q,id,j+1,N),N,n))nxt[mask|(1u<<j)]+=c;
          if(ti<T){uint w=(ti==0?sv2(pref,q,id,0,N):sv2(pref,q,id,d+ti-1,N));bool born=!nearFast(rank,id,d,mask,r,w,N,n);if((w==v&&born)||(w!=v&&!born))nxt[mask]+=c;}
        }
        for(uint z=0;z<4;z++)cur[z]=nxt[z];
      }
      ans+=cur[full];
    }
    return ans;
  }
  kernel void inverseRanks(device S2*x[[buffer(0)]],constant P2&p[[buffer(1)]],constant uint&np[[buffer(2)]],device uint*dead[[buffer(3)]],device const ushort*pref[[buffer(4)]],device ushort*rank[[buffer(5)]],uint id[[thread_position_in_grid]]){
    if(id>=np)return;for(uint v=0;v<p.N;v++)rank[id*MAXN2+v]=ushort(p.N);
    if(dead[0]||p.tau!=p.N+1)return;thread S2 q=x[id];
    for(uint k=0;k+1<p.N;k++)rank[id*MAXN2+uint(pref[id*MAXN2+k])]=ushort(k);
    for(uint v=0;v<p.N;v++)if(!has2(q,v)){rank[id*MAXN2+v]=ushort(p.N-1);break;}
  }
  kernel void inverseParticle(device S2*x[[buffer(0)]],constant P2&p[[buffer(1)]],constant uint&np[[buffer(2)]],device uint*dead[[buffer(3)]],device const ushort*pref[[buffer(4)]],device const ushort*rank[[buffer(5)]],device ulong*m2[[buffer(6)]],device ulong*m3[[buffer(7)]],uint id[[thread_position_in_grid]]){
    if(id>=np)return;if(dead[0]||p.tau!=p.N+1){m2[id]=0;m3[id]=0;return;}thread S2 q=x[id];m2[id]=invFast(pref,rank,q,id,p.N,p.n,2);m3[id]=invFast(pref,rank,q,id,p.N,p.n,3);
  }
  kernel void inverseReduce(device uint*a[[buffer(0)]],device const ulong*part2[[buffer(1)]],device const ulong*part3[[buffer(2)]],device ulong*m2[[buffer(3)]],device ulong*m3[[buffer(4)]],constant uint&np[[buffer(5)]],device uint*dead[[buffer(6)]],uint j[[thread_position_in_grid]]){
    if(dead[0]){m2[j]=0;m3[j]=0;return;}uint lo=256*j,hi=min(lo+256,np);ulong z2=0,z3=0;for(uint i=lo;i<hi;i++){z2+=ulong(a[i])*part2[i];z3+=ulong(a[i])*part3[i];}m2[j]=z2;m3[j]=z3;
  }
  """
// A second terminal-only library.  It consumes the already generated ordered
// prefix and cannot alter the SMC state, RNG, resampling, or inverse-fibre DP.
let geodesicMetal = """
  #include <metal_stdlib>
  using namespace metal;
  constant uint PREFIX_STRIDE=1024;
  struct P3{uint n,N,tau,r,lo,hi;};
  inline uint md3(uint a,uint b,uint n){return uint(abs(int(a%n)-int(b%n))+abs(int(a/n)-int(b/n)));}
  inline bool ed3(uint a,uint b,uint n){return md3(a,b,n)==1;}
  // One GPU thread owns one 256-particle group and writes a count for every
  // d=1,...,2n-1.  The prefix is terminal-ordered; d<=2n-1<=N-1, including
  // equality for n=2, so no synthetic final vertex is ever read.
  kernel void geodesicPrefixGroups(device const ushort*pref[[buffer(0)]],constant P3&p[[buffer(1)]],constant uint&np[[buffer(2)]],device const uint*dead[[buffer(3)]],device ulong*out[[buffer(4)]],uint j[[thread_position_in_grid]]){
    uint m=2u*p.n-1u,lo=256u*j,hi=min(lo+256u,np);
    if(dead[0]){for(uint d=0;d<m;d++)out[j*m+d]=0;return;}
    for(uint d=1;d<=m;d++){
      ulong total=0;
      for(uint id=lo;id<hi;id++){
        uint first=uint(pref[id*PREFIX_STRIDE]),last=first;bool good=true;
        for(uint k=1;k<d;k++){uint cur=uint(pref[id*PREFIX_STRIDE+k]);if(!ed3(last,cur,p.n))good=false;last=cur;}
        if(good&&md3(first,last,p.n)==d-1u)total++;
      }
      out[j*m+(d-1u)]=total;
    }
  }
  """
struct P { var n, N, tau, r, lo, hi: UInt32 }
struct R {
  let logZ, z, u, q, d, f, g, fOverZ, gOverZ, v, c, zm2, zm3, meanM2, meanM3: Double
  let m2Groups, m3Groups: [UInt64]
  let geodesicGroups: [UInt64]
  let terminalTotal: UInt64
  let dead: Bool
  let gpu, wall: Double
}
let dev = MTLCreateSystemDefaultDevice()!
let lib = try dev.makeLibrary(source: metal, options: nil)
let inverseLib = try dev.makeLibrary(source: inverseMetal, options: nil)
let geodesicLib = try dev.makeLibrary(source: geodesicMetal, options: nil)
let q = dev.makeCommandQueue()!
func pipe(_ s: String) -> MTLComputePipelineState {
  try! dev.makeComputePipelineState(function: lib.makeFunction(name: s)!)
}
let pc = pipe("count")
let pg = pipe("sum")
let pd = pipe("cdf")
let pp = pipe("prop")
let pm = pipe("move")
let pt = pipe("term")
let pir = try! dev.makeComputePipelineState(function: inverseLib.makeFunction(name: "inverseRanks")!)
let pip = try! dev.makeComputePipelineState(function: inverseLib.makeFunction(name: "inverseParticle")!)
let pig = try! dev.makeComputePipelineState(function: inverseLib.makeFunction(name: "inverseReduce")!)
let pgeo = try! dev.makeComputePipelineState(function: geodesicLib.makeFunction(name: "geodesicPrefixGroups")!)
var traceOn = false
var traceRows = [[String: Any]]()
func en(_ cb: MTLCommandBuffer, _ p: MTLComputePipelineState, _ bs: [(MTLBuffer, Int)], _ n: Int) {
  let e = cb.makeComputeCommandEncoder()!
  e.setComputePipelineState(p)
  for (b, i) in bs { e.setBuffer(b, offset: 0, index: i) }
  e.dispatchThreads(
    MTLSize(width: n, height: 1, depth: 1),
    threadsPerThreadgroup: MTLSize(width: 256, height: 1, depth: 1))
  e.endEncoding()
}
func run(_ n: Int, _ tau: Int, _ Np: Int, _ sd: UInt64, _ strat: Bool, _ moveL: Int) -> R {
  let N = n * n
  let G = (Np + 255) / 256
  var np = UInt32(Np)
  var gg = UInt32(G)
  var dead: UInt32 = 0
  var sm: UInt32 = strat ? 1 : 0
  var lm: UInt32 = UInt32(moveL)
  let nb = dev.makeBuffer(bytes: &np, length: 4, options: .storageModeShared)!
  let gbb = dev.makeBuffer(bytes: &gg, length: 4, options: .storageModeShared)!
  let db = dev.makeBuffer(bytes: &dead, length: 4, options: .storageModeShared)!
  let sb = dev.makeBuffer(bytes: &sm, length: 4, options: .storageModeShared)!
  let lb = dev.makeBuffer(bytes: &lm, length: 4, options: .storageModeShared)!
  let a = dev.makeBuffer(length: Np * 4, options: .storageModeShared)!
  let g = dev.makeBuffer(length: G * 8, options: .storageModeShared)!
  let c = dev.makeBuffer(length: G * 8, options: .storageModeShared)!
  let t = dev.makeBuffer(length: N * 8, options: .storageModeShared)!
  let mo = dev.makeBuffer(length: G * 8, options: .storageModeShared)!
  let inverseM2 = dev.makeBuffer(length: G * 8, options: .storageModeShared)!
  let inverseM3 = dev.makeBuffer(length: G * 8, options: .storageModeShared)!
  let geoM = 2 * n - 1
  // G*geoM <= 16384*63, so this terminal-only UInt64 buffer is at most
  // 8,257,536 bytes under the checked P<=2^22 bound.
  let geodesicGroups = dev.makeBuffer(length: G * geoM * 8, options: .storageModeShared)!
  // rank and per-particle inverse buffers belong solely to the terminal kernel.
  // At P=4096 this adds 8 MiB + 32 KiB + 32 KiB.
  let inverseRank = dev.makeBuffer(length: Np * 1024 * 2, options: .storageModeShared)!
  let inversePart2 = dev.makeBuffer(length: Np * 8, options: .storageModeShared)!
  let inversePart3 = dev.makeBuffer(length: Np * 8, options: .storageModeShared)!
  let rootF = dev.makeBuffer(length: G * 8, options: .storageModeShared)!
  let rootG = dev.makeBuffer(length: G * 8, options: .storageModeShared)!
  let pa = dev.makeBuffer(length: Np * 4, options: .storageModeShared)!
  var px = dev.makeBuffer(length: Np * 1024 * 2, options: .storageModeShared)!
  var py = dev.makeBuffer(length: Np * 1024 * 2, options: .storageModeShared)!
  let x = dev.makeBuffer(length: Np * 264, options: .storageModeShared)!
  let y = dev.makeBuffer(length: Np * 264, options: .storageModeShared)!
  memset(x.contents(), 0, Np * 264)
  memset(y.contents(), 0, Np * 264)
  memset(px.contents(), 0, Np * 1024 * 2)
  memset(py.contents(), 0, Np * 1024 * 2)
  var ps = [MTLBuffer]()
  var rs = [MTLBuffer]()
  for z in 1...N {
    var p = P(
      n: UInt32(n), N: UInt32(N), tau: UInt32(tau), r: UInt32(z),
      lo: UInt32(truncatingIfNeeded: sd), hi: UInt32(truncatingIfNeeded: sd >> 32))
    var r = UInt32(z - 1)
    ps.append(dev.makeBuffer(bytes: &p, length: 24, options: .storageModeShared)!)
    rs.append(dev.makeBuffer(bytes: &r, length: 4, options: .storageModeShared)!)
  }
  let st = Date()
  var gpuTotal = 0.0
  var src = x
  var dst = y
  for i in 0..<N {
    let cb = q.makeCommandBuffer()!
    en(cb, pc, [(src, 0), (a, 1), (ps[i], 2), (db, 3)], Np)
    en(cb, pg, [(a, 0), (g, 1), (nb, 2)], G)
    en(cb, pd, [(g, 0), (c, 1), (t, 2), (db, 3), (gbb, 4), (rs[i], 5)], 1)
    if i + 1 == N {
      en(cb, pt, [(src, 0), (a, 1), (mo, 2), (rootF, 3), (rootG, 4), (ps[i], 5), (nb, 6), (db, 7)], G)
      en(cb, pir, [(src, 0), (ps[i], 1), (nb, 2), (db, 3), (px, 4), (inverseRank, 5)], Np)
      en(cb, pip, [(src, 0), (ps[i], 1), (nb, 2), (db, 3), (px, 4), (inverseRank, 5), (inversePart2, 6), (inversePart3, 7)], Np)
      en(cb, pig, [(a, 0), (inversePart2, 1), (inversePart3, 2), (inverseM2, 3), (inverseM3, 4), (nb, 5), (db, 6)], G)
      en(cb, pgeo, [(px, 0), (ps[i], 1), (nb, 2), (db, 3), (geodesicGroups, 4)], G)
      cb.commit(); cb.waitUntilCompleted()
      precondition(cb.status == .completed && cb.error == nil)
      gpuTotal += cb.gpuEndTime-cb.gpuStartTime
      if traceOn {
        let aa=Array(UnsafeBufferPointer(start:a.contents().bindMemory(to:UInt32.self,capacity:Np),count:Np)).map(Int.init)
        let cc=Array(UnsafeBufferPointer(start:c.contents().bindMemory(to:UInt64.self,capacity:G),count:G)).map{String($0)}
        let words=(0..<Np).map { j -> [Int] in let p=px.contents().advanced(by:j*2048).bindMemory(to:UInt16.self,capacity:1024); return (0..<i).map{Int(p[$0])} }
        let bits=(0..<Np).map { j -> [Any] in let s=src.contents().advanced(by:j*264);let o=(0..<32).map{k in Int(s.load(fromByteOffset:k*4,as:UInt32.self))};let b=(0..<32).map{k in Int(s.load(fromByteOffset:128+k*4,as:UInt32.self))};return [o,b,Int(s.load(fromByteOffset:256,as:UInt32.self)),Int(s.load(fromByteOffset:260,as:UInt32.self))] }
        let mm=Array(UnsafeBufferPointer(start:mo.contents().bindMemory(to:UInt64.self,capacity:G),count:G)).map{String($0)}
        let ff=Array(UnsafeBufferPointer(start:rootF.contents().bindMemory(to:UInt64.self,capacity:G),count:G)).map{String($0)}
        let hh=Array(UnsafeBufferPointer(start:rootG.contents().bindMemory(to:Int64.self,capacity:G),count:G)).map{String($0)}
        let im2=Array(UnsafeBufferPointer(start:inverseM2.contents().bindMemory(to:UInt64.self,capacity:G),count:G)).map{String($0)}
        let im3=Array(UnsafeBufferPointer(start:inverseM3.contents().bindMemory(to:UInt64.self,capacity:G),count:G)).map{String($0)}
        let gg=geodesicGroups.contents().bindMemory(to:UInt64.self,capacity:G*geoM)
        let gp=(0..<G).map{j in (0..<geoM).map{k in String(gg[j*geoM+k])}}
        traceRows.append(["rank":i+1,"terminal":true,"allowed":aa,"cdf":cc,"prefix":words,"state":bits,"moment_groups":mm,"root_F_groups":ff,"root_G_groups":hh,"inverse_m2_groups":im2,"inverse_m3_groups":im3,"geodesic_prefix_groups":gp])
      }
    } else {
      en(
        cb, pp,
        [
          (src, 0), (a, 1), (c, 2), (t, 3), (dst, 4), (ps[i], 5), (nb, 6), (rs[i], 7), (db, 8),
          (sb, 9), (pa, 10), (px, 11), (py, 12),
        ], Np)
      if moveL > 0 { en(cb, pm, [(dst, 0), (py, 1), (ps[i], 2), (nb, 3), (lb, 4), (db, 5)], Np) }
      cb.commit(); cb.waitUntilCompleted()
      precondition(cb.status == .completed && cb.error == nil)
      gpuTotal += cb.gpuEndTime-cb.gpuStartTime
      if traceOn {
        let aa=Array(UnsafeBufferPointer(start:a.contents().bindMemory(to:UInt32.self,capacity:Np),count:Np)).map(Int.init)
        let cc=Array(UnsafeBufferPointer(start:c.contents().bindMemory(to:UInt64.self,capacity:G),count:G)).map{String($0)}
        let ppv=Array(UnsafeBufferPointer(start:pa.contents().bindMemory(to:UInt32.self,capacity:Np),count:Np)).map(Int.init)
        let words=(0..<Np).map { j -> [Int] in let p=py.contents().advanced(by:j*2048).bindMemory(to:UInt16.self,capacity:1024); return (0...i).map{Int(p[$0])} }
        let bits=(0..<Np).map { j -> [Any] in let s=dst.contents().advanced(by:j*264);let o=(0..<32).map{k in Int(s.load(fromByteOffset:k*4,as:UInt32.self))};let b=(0..<32).map{k in Int(s.load(fromByteOffset:128+k*4,as:UInt32.self))};return [o,b,Int(s.load(fromByteOffset:256,as:UInt32.self)),Int(s.load(fromByteOffset:260,as:UInt32.self))] }
        traceRows.append(["rank":i+1,"allowed":aa,"cdf":cc,"parents":ppv,"prefix":words,"state":bits])
      }
      swap(&src, &dst)
      swap(&px, &py)
    }
  }
  let ts = Array(
    UnsafeBufferPointer(start: t.contents().bindMemory(to: UInt64.self, capacity: N), count: N))
  if ts.contains(0) {
    return R(
      logZ: -.infinity, z: 0, u: 0, q: 0, d: 0, f: 0, g: 0, fOverZ: 0, gOverZ: 0, v: 0, c: 0, zm2: 0, zm3: 0, meanM2: 0, meanM3: 0, m2Groups: [], m3Groups: [], geodesicGroups: Array(repeating: 0, count: G * geoM), terminalTotal: 0, dead: true, gpu: gpuTotal,
      wall: -st.timeIntervalSinceNow)
  }
  var l = 0.0
  for (i, z) in ts.enumerated() { l += log(Double(z) / (Double(Np) * Double(N - i))) }
  let raw = UnsafeBufferPointer(
    start: mo.contents().bindMemory(to: UInt64.self, capacity: G), count: G
  ).reduce(0, &+)
  let rawF = UnsafeBufferPointer(
    start: rootF.contents().bindMemory(to: UInt64.self, capacity: G), count: G
  ).reduce(0, &+)
  let rawG = UnsafeBufferPointer(
    start: rootG.contents().bindMemory(to: Int64.self, capacity: G), count: G
  ).reduce(0, &+)
  let groupM2 = Array(UnsafeBufferPointer(
    start: inverseM2.contents().bindMemory(to: UInt64.self, capacity: G), count: G
  ))
  let groupM3 = Array(UnsafeBufferPointer(
    start: inverseM3.contents().bindMemory(to: UInt64.self, capacity: G), count: G
  ))
  let rawM2 = groupM2.reduce(0, &+)
  let rawM3 = groupM3.reduce(0, &+)
  let groupedGeodesic = Array(UnsafeBufferPointer(
    start: geodesicGroups.contents().bindMemory(to: UInt64.self, capacity: G * geoM), count: G * geoM
  ))
  let qv = Double(raw) / Double(ts[N - 1]) / Double(N)
  let z = exp(l)
  let u = z * Double(raw) / Double(ts[N - 1]) / Double(N)
  let d = z * Double(raw) / Double(ts[N - 1])
  let fOverZ = Double(rawF) / Double(ts[N - 1])
  let gOverZ = Double(rawG) / Double(ts[N - 1])
  let f = z * fOverZ
  let gValue = z * gOverZ
  return R(
    logZ: l, z: z, u: u, q: qv, d: d, f: f, g: gValue, fOverZ: fOverZ, gOverZ: gOverZ, v: fOverZ / 8, c: gOverZ / 4, zm2: z * Double(rawM2) / Double(ts[N - 1]), zm3: z * Double(rawM3) / Double(ts[N - 1]), meanM2: Double(rawM2) / Double(ts[N - 1]), meanM3: Double(rawM3) / Double(ts[N - 1]), m2Groups: groupM2, m3Groups: groupM3, geodesicGroups: groupedGeodesic, terminalTotal: ts[N - 1], dead: false, gpu: gpuTotal,
    wall: -st.timeIntervalSinceNow)
}
let v = CommandLine.arguments
let maybeMoveL = v.count >= 10 ? Int(v[9]) : 0
guard (v.count == 9 || v.count == 10 || v.count == 11), let n = Int(v[1]), let np = Int(v[2]), let Rn = Int(v[3]), let lo = Int(v[4]),
  let hi = Int(v[5]), let dom = UInt64(v[7], radix: 16), let moveL = maybeMoveL
else { fatalError("usage: bitset_smc n P R lo hi out domainhex strat|multi") }
guard v[8] == "strat" || v[8] == "multi" else { fatalError("method must be strat or multi") }
guard v.count < 11 || v[10] == "trace" else { fatalError("the optional final argument must be trace") }
let strat = v[8] == "strat"
traceOn = v.count == 11 && v[10] == "trace"
precondition(
  n >= 2 && n <= 32 && np > 0 && Rn > 0 && moveL >= 0 && lo >= 2 && hi <= n * n + 1 && lo <= hi
    // Every prefix and terminal-rank address is particle*1024+vertex in a
    // Metal uint.  This bound is checked before run() allocates any buffer.
    && np <= (1 << 22)
    && UInt64(np) <= UInt64.max / UInt64(np) / UInt64(n * n))
precondition(!traceOn || (Rn == 1 && lo == hi), "trace requires exactly one tau and one replicate")
let outURL = URL(fileURLWithPath: v[6])
FileManager.default.createFile(atPath: outURL.path, contents: nil)
let fh = try FileHandle(forWritingTo: outURL)
func o(_ s: String) { try! fh.write(contentsOf: (s + "\n").data(using: .utf8)!) }
o(
  "{\"kind\":\"metadata\",\"implementation\":\"r022-n32-terminal-geodesic-prefix\",\"base_source_sha256\":\"fc45ad181f6de297001c835a77d9be4620b944169269ed3cdb0f9c9a9a5a8629\",\"inverse_scope\":\"tau=N+1 only; canonical greedy inverse fibres d=2,3\",\"geodesic_scope\":\"terminal ordered prefix d=1,...,2n-1; every consecutive pair adjacent and endpoint Manhattan distance d-1\",\"n\":\(n),\"particles\":\(np),\"replicates\":\(Rn),\"tau_lo\":\(lo),\"tau_hi\":\(hi),\"domain\":\"\(String(dom,radix:16))\",\"method\":\"\(strat ? "strat" : "multi")\",\"move_steps\":\(moveL),\"moment_scale\":\"U_over_Z\",\"geodesic_integer_bound\":\"indicator in {0,1}; group<=256; raw<=P<=2^22; G*(2n-1)<=16384*63=1032192<2^32; bytes<=8257536\",\"inverse_integer_bound\":\"N<=1024: m2<=4092,m3<=4182024,a_i<=1; group<=256*m3=1070598144; host raw<=P*m3<=17540679991296<2^64 under P<=2^22\"}"
)
let T = Date()
for tau in lo...hi {
  for r in 0..<Rn {
    let sd = (0xA022_0000_0000_0000 &+ dom) &+ UInt64(n * 1_000_000 + tau * 10000 + r)
    let z: R = autoreleasepool { run(n, tau, np, sd, strat, moveL) }
    let l = z.dead ? "null" : String(z.logZ)
    let us = z.dead || tau > n * n ? "null" : String(z.u)
    let qs = z.dead || tau > n * n ? "null" : String(z.q)
    let rs = z.dead || tau > n * n
    let fs = rs ? "null" : String(z.f)
    let gs = rs ? "null" : String(z.g)
    let fzs = rs ? "null" : String(z.fOverZ)
    let gzs = rs ? "null" : String(z.gOverZ)
    let vs = rs ? "null" : String(z.v)
    let cs = rs ? "null" : String(z.c)
    let ds = rs ? "null" : String(z.d)
    let ids = rs ? "null" : String(z.f - 2 * z.g - 4 * z.d)
    let inverseActive = !z.dead && tau == n * n + 1
    let im2 = inverseActive ? String(z.zm2) : "null"
    let im3 = inverseActive ? String(z.zm3) : "null"
    let mm2 = inverseActive ? String(z.meanM2) : "null"
    let mm3 = inverseActive ? String(z.meanM3) : "null"
    let gm2 = inverseActive ? "[" + z.m2Groups.map(String.init).joined(separator: ",") + "]" : "null"
    let gm3 = inverseActive ? "[" + z.m3Groups.map(String.init).joined(separator: ",") + "]" : "null"
    let geoM = 2 * n - 1
    let geodesicActive = !z.dead
    let geoGroups = geodesicActive ? "[" + (0..<(z.geodesicGroups.count / geoM)).map { j in
      "[" + (0..<geoM).map { k in String(z.geodesicGroups[j * geoM + k]) }.joined(separator: ",") + "]"
    }.joined(separator: ",") + "]" : "null"
    let geoScaled = geodesicActive ? "[" + (0..<geoM).map { k in
      String(Double(stride(from: k, to: z.geodesicGroups.count, by: geoM).reduce(UInt64(0)) { $0 &+ z.geodesicGroups[$1] }) / Double(z.terminalTotal))
    }.joined(separator: ",") + "]" : "null"
    o(
      "{\"kind\":\"replicate\",\"tau\":\(tau),\"replicate\":\(r),\"seed_hex\":\"\(String(sd,radix:16))\",\"Z\":\(z.z),\"U\":\(us),\"logZ\":\(l),\"U_over_Z\":\(qs),\"scaled_U\":\(qs),\"D\":\(ds),\"F\":\(fs),\"G\":\(gs),\"F_over_Z\":\(fzs),\"G_over_Z\":\(gzs),\"V\":\(vs),\"C\":\(cs),\"F_minus_2G_minus_4D\":\(ids),\"inverse_m2_groups\":\(gm2),\"inverse_m3_groups\":\(gm3),\"Zmean_m2\":\(im2),\"Zmean_m3\":\(im3),\"mean_m2\":\(mm2),\"mean_m3\":\(mm3),\"scaled_mean_m2\":\(mm2),\"scaled_mean_m3\":\(mm3),\"geodesic_prefix_group_counts\":\(geoGroups),\"scaled_geodesic_prefix_p\":\(geoScaled),\"K2\":null,\"dead\":\(z.dead),\"gpu_seconds\":\(z.gpu),\"wall_seconds\":\(z.wall)}"
    )
    if traceOn {
      let obj:[String:Any] = ["implementation":"r022-n32-terminal-geodesic-prefix","base_source_sha256":"fc45ad181f6de297001c835a77d9be4620b944169269ed3cdb0f9c9a9a5a8629","method":strat ? "strat" : "multi","n":n,"tau":tau,"P":np,"seed_hex":String(sd,radix:16),"L":moveL,"rows":traceRows]
      let u=URL(fileURLWithPath:outURL.path+".trace.json")
      try JSONSerialization.data(withJSONObject:obj,options:[.prettyPrinted]).write(to:u)
      traceRows.removeAll()
    }
  }
}
o("{\"kind\":\"run_end\",\"wall_seconds\":\(-T.timeIntervalSinceNow)}")
