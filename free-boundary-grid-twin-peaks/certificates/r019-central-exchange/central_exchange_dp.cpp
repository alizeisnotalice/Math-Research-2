// Exact central-cut exchange audit for P_n square P_n, n <= 4.
// A word lists vertices from highest to lowest label.  For each source pair
// (p lower peak, q higher/global peak), q is therefore the forced first word.
// We count words having peak set {p,q} simultaneously in G and in the graph
// obtained by conjugating G with the transposition p <-> r.
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using Mask = uint32_t;
using Count = uint64_t;

struct Result {
  Count A2 = 0, same = 0, split = 0, same_oriented = 0, E = 0;
  int source_pairs = 0, positive_source_pairs = 0;
};

static std::vector<Mask> grid_neighbours(int n) {
  std::vector<Mask> nb(n * n);
  for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
    int v = y * n + x;
    for (auto [xx, yy] : {std::pair{x - 1, y}, {x + 1, y},
                           {x, y - 1}, {x, y + 1}})
      if (0 <= xx && xx < n && 0 <= yy && yy < n) nb[v] |= Mask(1) << (yy * n + xx);
  }
  return nb;
}

// t is a transposition, hence t^{-1}=t.  The returned graph is t^{-1} G t.
static std::vector<Mask> conjugate_swap(const std::vector<Mask>& g, int p, int r) {
  int N = static_cast<int>(g.size());
  std::vector<int> t(N);
  std::iota(t.begin(), t.end(), 0);
  std::swap(t[p], t[r]);
  std::vector<Mask> h(N);
  for (int v = 0; v < N; ++v) {
    Mask a = g[t[v]];
    while (a) {
      Mask b = a & -a; a -= b;
      h[v] |= Mask(1) << t[__builtin_ctz(b)];
    }
  }
  return h;
}

// q is the unique global maximum, hence first.  All subsequent vertices other
// than p must have an earlier neighbour in both graphs; p must have none in
// either graph.  This is exactly simultaneous peak set {p,q}.
static Count simultaneous_words(const std::vector<Mask>& g, const std::vector<Mask>& h,
                                 int p, int q) {
  const int N = static_cast<int>(g.size());
  const Mask full = (Mask(1) << N) - 1;
  std::vector<Count> dp(size_t(1) << N);
  dp[Mask(1) << q] = 1;
  for (Mask S = 0; S <= full; ++S) {
    Count ways = dp[S];
    if (!ways) continue;
    for (int v = 0; v < N; ++v) if (!(S & (Mask(1) << v))) {
      if (v == q) continue;
      bool g_birth = !(g[v] & S), h_birth = !(h[v] & S);
      bool ok = (v == p) ? (g_birth && h_birth) : (!g_birth && !h_birth);
      if (ok) dp[S | (Mask(1) << v)] += ways;
    }
  }
  return dp[full];
}

static bool is_peak_word(const std::vector<int>& word, const std::vector<Mask>& g,
                         int p, int q) {
  Mask S = 0, peaks = 0;
  for (int v : word) {
    if (!(g[v] & S)) peaks |= Mask(1) << v;
    S |= Mask(1) << v;
  }
  return peaks == ((Mask(1) << p) | (Mask(1) << q));
}

// Independent calibration for n=2,3: direct enumeration with q first.
static Count brute_simultaneous(const std::vector<Mask>& g, const std::vector<Mask>& h,
                                int p, int q) {
  const int N = static_cast<int>(g.size());
  std::vector<int> rest;
  for (int v = 0; v < N; ++v) if (v != q) rest.push_back(v);
  Count ans = 0;
  do {
    std::vector<int> word{q}; word.insert(word.end(), rest.begin(), rest.end());
    if (is_peak_word(word, g, p, q) && is_peak_word(word, h, p, q)) ++ans;
  } while (std::next_permutation(rest.begin(), rest.end()));
  return ans;
}

static int calibration(int n) {
  const auto g = grid_neighbours(n);
  const int N = n * n;
  int tested = 0;
  for (int p = 0; p < N; ++p) for (int q = 0; q < N; ++q) {
    if (p == q || (g[p] & (Mask(1) << q))) continue;
    for (int r = 0; r < N; ++r) {
      if (r == p || r == q) continue;
      auto h = conjugate_swap(g, p, r);
      Count a = simultaneous_words(g, h, p, q);
      Count b = brute_simultaneous(g, h, p, q);
      assert(a == b);
      ++tested;
    }
  }
  return tested;
}

static Result central_audit(int n) {
  assert(n % 2 == 0 && n <= 4);
  const auto g = grid_neighbours(n);
  const int N = n * n, m = n / 2;
  Result z;
  for (int p = 0; p < N; ++p) for (int q = p + 1; q < N; ++q) {
    if (g[p] & (Mask(1) << q)) continue;
    Count pair = simultaneous_words(g, g, p, q) + simultaneous_words(g, g, q, p);
    z.A2 += pair;
    bool same_side = ((p % n) < m) == ((q % n) < m);
    if (same_side) z.same += pair; else z.split += pair;
  }
  // Ordered (p,q): p is the lower peak and q is the global maximum.
  for (int p = 0; p < N; ++p) for (int q = 0; q < N; ++q) {
    if (p == q || (g[p] & (Mask(1) << q))) continue;
    bool same_side = ((p % n) < m) == ((q % n) < m);
    if (!same_side) continue;
    Count original_ordered = simultaneous_words(g, g, p, q);
    z.same_oriented += original_ordered;
    int x = p % n, y = p / n;
    int r = y * n + (n - 1 - x);
    assert(r != q);
    Count e = simultaneous_words(g, conjugate_swap(g, p, r), p, q);
    assert(e <= original_ordered);  // simultaneous words form a literal subset.
    z.E += e;
    ++z.source_pairs;
    if (e) ++z.positive_source_pairs;
  }
  assert(z.same_oriented == z.same);
  assert(z.A2 == z.same + z.split);
  assert(z.E <= z.split);  // required by the intended injection.
  return z;
}

int main() {
  int c2 = calibration(2), c3 = calibration(3);
  Result n4 = central_audit(4);
  // Known exact A_2 values provide a separate checksum for the ordinary DP.
  assert(n4.A2 == 779299174080ULL);
  std::cout << "{\"status\":\"passed\",\"calibration\":{\"n2_triples\":" << c2
            << ",\"n3_triples\":" << c3 << "},\"n4\":{"
            << "\"A2\":\"" << n4.A2 << "\",\"same\":\"" << n4.same
            << "\",\"split\":\"" << n4.split << "\",\"same_oriented\":\""
            << n4.same_oriented << "\",\"E\":\"" << n4.E
            << "\",\"same_source_pairs\":" << n4.source_pairs
            << ",\"positive_source_pairs\":" << n4.positive_source_pairs << "}}\n";
}
