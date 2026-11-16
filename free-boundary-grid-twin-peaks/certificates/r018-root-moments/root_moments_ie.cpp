// Exact K=2 root moments through the project's independent-set C(I) recurrence.
// The upstream source is included read-only; its main is renamed, and this
// file independently reruns the recurrence while adding signed moment sums.
#define main upstream_independent_peak_dp_main
#include "../../computation/independent_peak_dp.cpp"
#undef main

static void addsigned(U256& plus, U256& minus, const U192& v, int64_t coeff) {
  if (coeff >= 0) addscaled(plus, v, uint64_t(coeff));
  else addscaled(minus, v, uint64_t(-coeff));
}

int main(int argc, char** argv) {
  if (argc != 2) return 2;
  n = std::atoi(argv[1]); if (n < 2 || n > 6) return 2;
  N = n*n; alpha = (N+1)/2;
  valid.clear();
  for (uint64_t r=0; r<(1ULL<<n); ++r) if (!(r&(r<<1))) valid.push_back(r);
  bags.assign(alpha+1, {}); enumerate(0, 0, 0);
  size_t total=0; for (auto& b:bags) total += b.size();
  std::vector<uint64_t> sets; sets.reserve(total);
  std::vector<size_t> start(alpha+2);
  for (int k=0;k<=alpha;k++) {
    start[k]=sets.size(); sets.insert(sets.end(), bags[k].begin(), bags[k].end());
    std::vector<uint64_t>().swap(bags[k]);
  }
  start[alpha+1]=sets.size();
  std::unordered_map<uint64_t,uint32_t> index; index.max_load_factor(.7);
  index.reserve(total*2);
  for (uint32_t i=0;i<sets.size();i++) index.emplace(sets[i],i);
  std::vector<U192> C(total); C[0].x[0]=1;
  for (int i=2;i<=N;i++) C[0]=mulsmall192(C[0],i);
  for (int k=1;k<=alpha;k++) for(size_t z=start[k];z<start[k+1];z++) {
    uint64_t s=sets[z], q=s; U192 sum;
    while(q) { uint64_t bit=q&-q; q-=bit; sum=add192(sum,C[index.at(s^bit)]); }
    C[z]=divsmall192(sum,k+__builtin_popcountll(neighbors(s)));
  }

  U256 a2p,a2m,firstp,firstm,dotp,dotm;
  std::vector<U256> margp(N), margm(N);
  for (size_t z=0;z<sets.size();z++) {
    const U192& cval=C[z]; if (iszero(cval)) continue;
    uint64_t s=sets[z]; int k=__builtin_popcountll(s); if(k<2) continue;
    const bool ie_negative=(k&1); // (-1)^(k-2)
    addsigned(a2p,a2m,cval,ie_negative ? -int64_t(choose(k,2)) : int64_t(choose(k,2)));
    int64_t sx=0, sy=0, normsum=0; std::array<int,18> vs{}; int nv=0;
    for(uint64_t q=s;q;) {
      uint64_t bit=q&-q; q-=bit; int v=__builtin_ctzll(bit); vs[nv++]=v;
      int64_t tx=2*(v%n)-(n-1), ty=2*(v/n)-(n-1);
      sx+=tx; sy+=ty; normsum+=tx*tx+ty*ty;
      addsigned(margp[v],margm[v],cval,ie_negative ? -int64_t(k-1) : int64_t(k-1));
    }
    int64_t first=(k-1)*normsum;
    int64_t dot=(sx*sx+sy*sy-normsum)/2;
    assert((sx*sx+sy*sy-normsum)%2==0);
    addsigned(firstp,firstm,cval,ie_negative ? -first : first);
    addsigned(dotp,dotm,cval,ie_negative ? -dot : dot);
  }
  std::cout << "{\"n\":"<<n<<",\"N\":"<<N<<",\"independent_sets\":"<<total
    <<",\"A2\":\""<<signedstr(a2p,a2m)<<"\",\"firstsum\":\""<<signedstr(firstp,firstm)
    <<"\",\"secondsum\":\""<<signedstr(dotp,dotm)<<"\",\"root_marginals\":[";
  for(int v=0;v<N;v++) { if(v) std::cout<<','; std::cout<<'\"'<<signedstr(margp[v],margm[v])<<'\"'; }
  std::cout << "]}\n";
}
