// Exact ordered first-merge sector ledger on P_n, n=3 or 4 only.
//
// A vertex v is born only adjacent to the occupied set.  h(S) counts
// completions from S.  a[S,r] counts admissible histories of S whose first
// (low-root) vertex is r; it is obtained by deleting the final vertex.  At a
// first merge S=A union B union {w}, the low root belongs to exactly one of
// A,B.  We enumerate both directions explicitly and record b=|low branch|.
// All output accumulators are exact unsigned __int128 values.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using U128 = unsigned __int128;

static std::string dec(U128 x) {
  if (x == 0) return "0";
  std::string out;
  while (x) { out.push_back(char('0' + x % 10)); x /= 10; }
  std::reverse(out.begin(), out.end());
  return out;
}

static uint64_t choose64(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  uint64_t z = 1;
  for (int i = 1; i <= k; ++i) z = z * uint64_t(n - k + i) / uint64_t(i);
  return z;
}

static uint32_t component(uint32_t seed, uint32_t allowed,
                          const std::vector<uint32_t>& adj) {
  uint32_t seen = seed, frontier = seed;
  while (frontier) {
    const uint32_t bit = frontier & -frontier;
    frontier ^= bit;
    const int v = __builtin_ctz(bit);
    const uint32_t add = adj[v] & allowed & ~seen;
    seen |= add;
    frontier |= add;
  }
  return seen;
}

static int coordx(int v, int n) { return v % n; }
static int coordy(int v, int n) { return v / n; }

struct Cell { U128 count = 0, d2sum = 0, rho_gt3 = 0; };

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: ordered_sector_dp N OUTPUT.json\n";
    return 2;
  }
  const int n = std::stoi(argv[1]);
  if (n < 3 || n > 4) {
    std::cerr << "hard guard: only n=3,4 are authorized\n";
    return 2;
  }
  const int N = n * n;
  const uint32_t cap = uint32_t(1) << N;
  const uint32_t full = cap - 1;
  std::vector<uint32_t> adj(N, 0);
  for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
    const int v = y * n + x;
    for (const auto& d : std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}) {
      const int X = x + d[0], Y = y + d[1];
      if (0 <= X && X < n && 0 <= Y && Y < n) adj[v] |= uint32_t(1) << (Y*n+X);
    }
  }

  // 'connected' is explicit to make the a/h recurrence auditably scoped.
  std::vector<uint8_t> connected(cap, 0);
  connected[0] = 0;
  for (uint32_t S = 1; S < cap; ++S)
    connected[S] = (component(S & -S, S, adj) == S);

  // Rooted history table.  uint64_t is sufficient for every individual n=4
  // entry (the aggregate is nevertheless U128); this is 2^16*16*8 = 8 MiB.
  std::vector<uint64_t> a(size_t(cap) * N, 0);
  auto A = [&](uint32_t S, int r) -> uint64_t& { return a[size_t(S)*N + r]; };
  for (int r = 0; r < N; ++r) A(uint32_t(1) << r, r) = 1;
  for (int k = 2; k <= N; ++k) {
    for (uint32_t S = 1; S < cap; ++S) {
      if (__builtin_popcount(S) != k || !connected[S]) continue;
      for (uint32_t qr = S; qr; qr &= qr-1) {
        const int r = __builtin_ctz(qr);
        uint64_t z = 0;
        for (uint32_t qv = S ^ (uint32_t(1)<<r); qv; qv &= qv-1) {
          const uint32_t bit = qv & -qv;
          const uint32_t T = S ^ bit;
          const int v = __builtin_ctz(bit);
          if (adj[v] & T) z += A(T, r);
        }
        A(S, r) = z;
      }
    }
  }

  // Completion counts h(S), descending in |S|.
  std::vector<uint64_t> h(cap, 0);
  h[full] = 1;
  for (int k = N-1; k >= 1; --k) for (uint32_t S = 1; S < cap; ++S) {
    if (__builtin_popcount(S) != k || !connected[S]) continue;
    uint64_t z = 0;
    for (uint32_t q = full ^ S; q; q &= q-1) {
      const uint32_t bit = q & -q;
      if (adj[__builtin_ctz(bit)] & S) z += h[S | bit];
    }
    h[S] = z;
  }

  // Index q=a+b is the merger rank minus one; b is the *ordered* low branch.
  std::vector<Cell> by_qb(size_t(N+1) * (N+1));
  auto cell = [&](int q, int b) -> Cell& { return by_qb[size_t(q)*(N+1)+b]; };
  std::vector<U128> manhattan(2*(n-1)+1, 0);
  std::vector<U128> euclidean_squared(2*(n-1)*(n-1)+1, 0);
  U128 A1 = 0, A2 = 0, D2 = 0, rho_gt3 = 0;
  for (int r = 0; r < N; ++r) A1 += A(full, r);
  uint64_t eligible_states = 0, merge_cuts = 0;

  for (int s = 3; s <= N; ++s) for (uint32_t S = 1; S < cap; ++S) {
    if (__builtin_popcount(S) != s || !connected[S] || h[S] == 0) continue;
    ++eligible_states;
    for (uint32_t qw = S; qw; qw &= qw-1) {
      const uint32_t wbit = qw & -qw;
      const int w = __builtin_ctz(wbit);
      const uint32_t T = S ^ wbit;
      const uint32_t B0 = T & -T;
      const uint32_t B = component(B0, T, adj);
      const uint32_t C = T ^ B;
      if (C == 0 || component(C & -C, C, adj) != C) continue;
      if (!(adj[w] & B) || !(adj[w] & C)) continue;
      ++merge_cuts;
      const int bB = __builtin_popcount(B), bC = __builtin_popcount(C);
      const int qrank = s - 1;
      // Each component is selected as the low branch once.  Thus this is not
      // min(branch size), nor an unordered pair ledger.
      for (int direction = 0; direction < 2; ++direction) {
        const uint32_t L = direction == 0 ? B : C;
        const uint32_t H = direction == 0 ? C : B;
        const int low_size = direction == 0 ? bB : bC;
        const uint64_t interleave = choose64(qrank - 1, low_size);
        for (uint32_t qr = L; qr; qr &= qr-1) {
          const int r = __builtin_ctz(qr);
          const uint64_t ar = A(L, r);
          if (!ar) continue;
          for (uint32_t qp = H; qp; qp &= qp-1) {
            const int p = __builtin_ctz(qp);
            const uint64_t ap = A(H, p);
            if (!ap) continue;
            const U128 weight = U128(h[S]) * interleave * ar * ap;
            const int dx = coordx(r,n) - coordx(p,n);
            const int dy = coordy(r,n) - coordy(p,n);
            const int d2 = dx*dx + dy*dy;
            const int d1 = std::abs(dx) + std::abs(dy);
            Cell& z = cell(qrank, low_size);
            z.count += weight; z.d2sum += weight * d2;
            A2 += weight; D2 += weight * d2;
            manhattan[d1] += weight; euclidean_squared[d2] += weight;
            if (d1 > 3) { rho_gt3 += weight; z.rho_gt3 += weight; }
          }
        }
      }
    }
  }

  std::ofstream out(argv[2]);
  if (!out) { std::cerr << "cannot open output\n"; return 3; }
  out << "{\n"
      << "  \"schema\": \"ordered-first-merge-sector-v1\",\n"
      << "  \"n\": " << n << ",\n"
      << "  \"workers\": 1,\n"
      << "  \"hard_guard\": \"3 <= n <= 4\",\n"
      << "  \"rooted_table_bytes\": " << uint64_t(a.size()*sizeof(uint64_t)) << ",\n"
      << "  \"A1\": \"" << dec(A1) << "\",\n"
      << "  \"A2\": \"" << dec(A2) << "\",\n"
      << "  \"D_euclidean_squared_sum\": \"" << dec(D2) << "\",\n"
      << "  \"rho_manhattan_gt_3_count\": \"" << dec(rho_gt3) << "\",\n"
      << "  \"eligible_connected_states\": " << eligible_states << ",\n"
      << "  \"merge_cuts\": " << merge_cuts << ",\n"
      << "  \"ordered_q_lowb\": [\n";
  bool first = true;
  for (int q = 2; q < N; ++q) for (int b = 1; b < q; ++b) {
    const Cell& z = cell(q,b); if (z.count == 0) continue;
    if (!first) out << ",\n"; first = false;
    out << "    {\"q\": " << q << ", \"low_b\": " << b
        << ", \"count\": \"" << dec(z.count) << "\", \"D_euclidean_squared_sum\": \""
        << dec(z.d2sum) << "\", \"rho_manhattan_gt_3_count\": \"" << dec(z.rho_gt3) << "\"}";
  }
  out << "\n  ],\n  \"manhattan_distance_histogram\": [";
  for (size_t d = 0; d < manhattan.size(); ++d) {
    if (d) out << ", ";
    out << "{\"distance\": " << d << ", \"count\": \"" << dec(manhattan[d]) << "\"}";
  }
  out << "],\n  \"euclidean_squared_distance_histogram\": [";
  bool firste = true;
  for (size_t e = 0; e < euclidean_squared.size(); ++e) if (euclidean_squared[e]) {
    if (!firste) out << ", "; firste = false;
    out << "{\"squared_distance\": " << e << ", \"count\": \"" << dec(euclidean_squared[e]) << "\"}";
  }
  out << "]\n}\n";
  std::cout << "n=" << n << " A1=" << dec(A1) << " A2=" << dec(A2)
            << " D2=" << dec(D2) << " rho>3=" << dec(rho_gt3) << "\n";
  return 0;
}
