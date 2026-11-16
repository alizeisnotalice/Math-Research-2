// Exact ordered-domino promotion accounting through the upstream C(I) DP.
#define main upstream_independent_peak_dp_main
#include "../../computation/independent_peak_dp.cpp"
#undef main

static void addsigned(U256& plus, U256& minus, const U192& value, int64_t coeff) {
  if (coeff >= 0) addscaled(plus, value, uint64_t(coeff));
  else addscaled(minus, value, uint64_t(-coeff));
}

static int nbor_edges = 0;
static std::array<uint64_t, 36> closedv{};
static int ordered_domino_cover(uint64_t P) {
  int ans = 0;
  for (int a=0; a<N; ++a) {
    uint64_t ns = neighbors(1ULL << a);
    while (ns) {
      uint64_t bit = ns & -ns; ns -= bit;
      int b = __builtin_ctzll(bit);
      if (!(P & ~(closedv[a] | closedv[b]))) ++ans;
    }
  }
  return ans;
}
static int common_neighbour_count(uint64_t P) {
  uint64_t ans = (1ULL<<N)-1, q=P;
  while (q) { uint64_t bit=q&-q; q-=bit; ans &= neighbors(bit); }
  return __builtin_popcountll(ans);
}
static bool has_distance_three_pair(uint64_t P) {
  for (uint64_t q=P; q; q&=q-1) {
    int a=__builtin_ctzll(q);
    uint64_t tail=P & ~((1ULL<<(a+1))-1);
    while (tail) { uint64_t bit=tail&-tail; tail-=bit; int b=__builtin_ctzll(bit);
      if (std::abs(a%n-b%n)+std::abs(a/n-b/n)==3) return true;
    }
  }
  return false;
}

int main(int argc, char** argv) {
  if (argc != 2) return 2;
  n = std::atoi(argv[1]); if (n < 2 || n > 5) return 2;
  N = n*n; alpha = (N+1)/2;
  valid.clear();
  for (uint64_t r=0; r<(1ULL<<n); ++r) if (!(r&(r<<1))) valid.push_back(r);
  bags.assign(alpha+1, {}); enumerate(0, 0, 0);
  size_t total = 0; for (auto& bag : bags) total += bag.size();
  std::vector<uint64_t> sets; sets.reserve(total);
  std::vector<size_t> start(alpha+2);
  for (int k=0; k<=alpha; ++k) {
    start[k] = sets.size(); sets.insert(sets.end(), bags[k].begin(), bags[k].end());
    std::vector<uint64_t>().swap(bags[k]);
  }
  start[alpha+1] = sets.size();
  std::unordered_map<uint64_t,uint32_t> index; index.max_load_factor(.7); index.reserve(total*2);
  for (uint32_t i=0; i<sets.size(); ++i) index.emplace(sets[i], i);
  std::vector<U192> C(total); C[0].x[0] = 1;
  for (int i=2; i<=N; ++i) C[0] = mulsmall192(C[0], i);
  for (int k=1; k<=alpha; ++k) for (size_t z=start[k]; z<start[k+1]; ++z) {
    uint64_t s=sets[z], q=s; U192 sum;
    while (q) { uint64_t bit=q&-q; q-=bit; sum=add192(sum, C[index.at(s^bit)]); }
    C[z] = divsmall192(sum, k+__builtin_popcountll(neighbors(s)));
  }
  for (int v=0; v<N; ++v) closedv[v] = neighbors(1ULL<<v) | (1ULL<<v);

  std::vector<U256> wp(alpha+1), wm(alpha+1);
  std::vector<U256> nonstarp(alpha+1), nonstarm(alpha+1);
  std::vector<U256> d2p(2*(n-1)+1), d2m(2*(n-1)+1);
  U256 lhsp, lhsm, rhsp, rhsm, starp, starm, star2p, star2m;
  for (size_t z=0; z<sets.size(); ++z) {
    uint64_t I=sets[z]; int k=__builtin_popcountll(I); const U192& c=C[z];
    if (iszero(c)) continue;
    // LHS: (N)_2 times the ordinary exact-one-peak IE coefficient.
    int64_t lcoeff = int64_t(N) * (N-1) * k;
    if ((k-1)&1) lcoeff = -lcoeff;
    addsigned(lhsp, lhsm, c, lcoeff);

    // Möbius transform every peak-subset weight e_ord(P).
    for (uint64_t P=I;; P=(P-1)&I) {
      int r=__builtin_popcountll(P);
      if (r) {
        int e=ordered_domino_cover(P);
        if (e) {
          int64_t coeff=e;
          if ((k-r)&1) coeff=-coeff;
          addsigned(wp[r], wm[r], c, coeff);
          addsigned(rhsp, rhsm, c, coeff);
          if (common_neighbour_count(P)) {
            addsigned(starp, starm, c, coeff);
            if (r>=2) addsigned(star2p, star2m, c, coeff);
          }
          else {
            assert(has_distance_three_pair(P));
            addsigned(nonstarp[r], nonstarm[r], c, coeff);
          }
          if (r==2) {
            int a=__builtin_ctzll(P), b=__builtin_ctzll(P & ~(1ULL<<a));
            int d=std::abs(a%n-b%n)+std::abs(a/n-b/n);
            addsigned(d2p[d], d2m[d], c, coeff);
          }
        }
      }
      if (!P) break;
    }
  }
  assert(signedstr(lhsp,lhsm) == signedstr(rhsp,rhsm));
  std::cout << "{\"status\":\"passed\",\"n\":"<<n<<",\"N_falling_2_A1\":\""
            << signedstr(lhsp,lhsm) << "\",\"RHS\":\"" << signedstr(rhsp,rhsm)
            << "\",\"weighted_by_peak_count\":{";
  for (int r=1;r<=alpha;++r) { if(r>1)std::cout<<','; std::cout<<'"'<<r<<"\":\""<<signedstr(wp[r],wm[r])<<'"'; }
  std::cout << "},\"star_common_neighbour_weight\":\""<<signedstr(starp,starm)
            << "\",\"star_K_ge_2_weight\":\""<<signedstr(star2p,star2m)
            << "\",\"nonstar_weighted_by_peak_count\":{";
  for (int r=1;r<=alpha;++r) { if(r>1)std::cout<<','; std::cout<<'"'<<r<<"\":\""<<signedstr(nonstarp[r],nonstarm[r])<<'"'; }
  std::cout << "},\"K2_weighted_by_manhattan_distance\":{";
  bool first=true; for (int d=0;d<=2*(n-1);++d) if (cmp(d2p[d],d2m[d])) {
    if(!first)std::cout<<','; first=false; std::cout<<'"'<<d<<"\":\""<<signedstr(d2p[d],d2m[d])<<'"';
  }
  std::cout << "}}\n";
}
