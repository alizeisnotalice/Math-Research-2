// Exact K=2 joint count by peak distance and lower-peak descending birth rank.
// The high peak q is the global maximum and is forced first. Before the lower
// peak p is inserted, every added nonroot must touch the earlier set; p is
// inserted only with no earlier neighbour, at tau=|S|+1.  After p, the same
// attachment rule applies to every remaining vertex.
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using U64 = std::uint64_t;

static void add_checked(U64& a, U64 b) {
  assert(UINT64_MAX - a >= b);
  a += b;
}
static void add_product(U64& a, U64 x, U64 y) {
  __uint128_t z = (__uint128_t)x * y + a;
  assert(z <= UINT64_MAX);
  a = (U64)z;
}

int main(int argc, char** argv) {
  if (argc != 2) return 2;
  const int n = std::atoi(argv[1]);
  if (n < 2 || n > 4) return 2;
  const int N = n*n, ALL = (1 << N) - 1;
  std::array<int, 16> nbr{};
  for (int y=0; y<n; ++y) for (int x=0; x<n; ++x) {
    const int v=y*n+x;
    for (auto [xx, yy] : {std::pair{x-1,y}, {x+1,y}, {x,y-1}, {x,y+1}})
      if (0<=xx && xx<n && 0<=yy && yy<n) nbr[v] |= 1 << (yy*n+xx);
  }
  std::array<std::array<U64,17>, 17> joint{};
  std::vector<U64> suffix(1 << N), prefix(1 << N);
  for (int p=0; p<N; ++p) {
    const int pb = 1 << p;
    // suffix[S] is the number of valid completions after p has already been
    // added. Descending numerical-mask order is valid because each move sets
    // one bit and hence increases the mask.
    suffix[ALL] = 1;
    for (int S=ALL-1; S>=0; --S) {
      if (!(S & pb)) { suffix[S]=0; continue; }
      U64 total=0;
      for (int v=0; v<N; ++v) if (!(S&(1<<v)) && (S&nbr[v]))
        add_checked(total, suffix[S|(1<<v)]);
      suffix[S]=total;
    }
    for (int q=0; q<N; ++q) if (q != p) {
      const int qb=1<<q;
      std::fill(prefix.begin(), prefix.end(), 0);
      prefix[qb]=1;
      for (int S=0; S<=ALL; ++S) {
        const U64 ways=prefix[S];
        if (!ways || (S&pb)) continue;
        // Insert the lower peak p now iff it has no earlier neighbour.
        if (!(S&nbr[p])) {
          const int tau=__builtin_popcount((unsigned)S)+1;
          const int d=std::abs(p%n-q%n)+std::abs(p/n-q/n);
          add_product(joint[d][tau], ways, suffix[S|pb]);
        }
        // A nonroot cannot be an additional birth.
        for (int v=0; v<N; ++v) if (v!=p && !(S&(1<<v)) && (S&nbr[v]))
          add_checked(prefix[S|(1<<v)], ways);
      }
    }
  }
  U64 total=0;
  std::array<U64,17> dh{}, th{};
  for(int d=0;d<=2*(n-1);++d) for(int t=2;t<=N;++t) {
    add_checked(total,joint[d][t]); add_checked(dh[d],joint[d][t]); add_checked(th[t],joint[d][t]);
  }
  std::cout << "{\"n\":"<<n<<",\"N\":"<<N
            <<",\"method\":\"fixed ordered roots, prefix/suffix subset DP\""
            <<",\"joint_distance_birth\":{";
  bool firstd=true;
  for(int d=1;d<=2*(n-1);++d) if(dh[d]) {if(!firstd)std::cout<<',';firstd=false;std::cout<<'\"'<<d<<"\":{";bool firstt=true;for(int t=2;t<=N;++t)if(joint[d][t]){if(!firstt)std::cout<<',';firstt=false;std::cout<<'\"'<<t<<"\":\""<<joint[d][t]<<'\"';}std::cout<<'}';}
  std::cout << "},\"distance_counts\":{";bool comma=false;for(int d=1;d<=2*(n-1);++d)if(dh[d]){if(comma)std::cout<<',';comma=true;std::cout<<'\"'<<d<<"\":\""<<dh[d]<<'\"';}
  std::cout << "},\"birth_counts\":{";comma=false;for(int t=2;t<=N;++t)if(th[t]){if(comma)std::cout<<',';comma=true;std::cout<<'\"'<<t<<"\":\""<<th[t]<<'\"';}
  std::cout << "},\"A2\":\""<<total<<"\"}\n";
}
