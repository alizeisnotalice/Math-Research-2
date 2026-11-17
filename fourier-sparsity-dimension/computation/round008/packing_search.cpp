// Bounded exact search for disjoint affine codimension-3 flats in F_2^n.
// A pilot uses a depth-3 decision-tree partition with seven independent
// query normals, punctures 1--4 leaves, then makes exact bitset-valid
// add/remove/replace moves.  A hit has collective normal rank >= target.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using std::array;
using std::uint16_t;
using std::uint64_t;
using std::vector;

constexpr int MAXWORDS = 8; // n <= 9, so 2^n <= 512

struct Flat {
  array<uint16_t, 3> normal{};
  array<unsigned char, 3> rhs{};
  array<uint64_t, MAXWORDS> point{};
};

static int parity(uint16_t x) { return __builtin_parity((unsigned)x); }

static int gf2rank(vector<uint16_t> a, int n) {
  int r = 0;
  for (int b = n - 1; b >= 0; --b) {
    int p = r;
    while (p < (int)a.size() && !(a[p] & (uint16_t(1) << b))) ++p;
    if (p == (int)a.size()) continue;
    std::swap(a[r], a[p]);
    for (int i = 0; i < (int)a.size(); ++i)
      if (i != r && (a[i] & (uint16_t(1) << b))) a[i] ^= a[r];
    ++r;
  }
  return r;
}

static int packing_rank(const vector<Flat>& p, int n) {
  vector<uint16_t> rows;
  for (const Flat& f : p) rows.insert(rows.end(), f.normal.begin(), f.normal.end());
  return gf2rank(rows, n);
}

static Flat make_flat(array<uint16_t, 3> a, array<unsigned char, 3> b, int n) {
  if (gf2rank({a[0], a[1], a[2]}, n) != 3) throw std::runtime_error("dependent normal rows");
  Flat f; f.normal = a; f.rhs = b;
  const int N = 1 << n;
  for (int x = 0; x < N; ++x) {
    if (parity(a[0] & x) == b[0] && parity(a[1] & x) == b[1] && parity(a[2] & x) == b[2])
      f.point[x >> 6] |= uint64_t(1) << (x & 63);
  }
  return f;
}

static bool disjoint(const Flat& f, const vector<Flat>& p, int skip = -1) {
  for (int i = 0; i < (int)p.size(); ++i) if (i != skip)
    for (int w = 0; w < MAXWORDS; ++w) if (f.point[w] & p[i].point[w]) return false;
  return true;
}

static Flat random_flat(std::mt19937_64& rng, int n) {
  const uint16_t lim = uint16_t(1u << n);
  array<uint16_t,3> a{};
  do {
    for (auto& x : a) x = uint16_t(1 + rng() % (lim - 1));
  } while (gf2rank({a[0], a[1], a[2]}, n) != 3);
  return make_flat(a, {static_cast<unsigned char>(rng() & 1), static_cast<unsigned char>(rng() & 1), static_cast<unsigned char>(rng() & 1)}, n);
}

static vector<uint16_t> independent_queries(std::mt19937_64& rng, int n) {
  vector<uint16_t> q;
  const uint16_t lim = uint16_t(1u << n);
  while ((int)q.size() < 7) {
    uint16_t x = uint16_t(1 + rng() % (lim - 1));
    vector<uint16_t> t = q; t.push_back(x);
    if (gf2rank(t, n) > (int)q.size()) q.push_back(x);
  }
  return q;
}

// Leaves of a depth-three binary decision tree, giving an exact eight-flat partition.
static vector<Flat> decision_tree_partition(std::mt19937_64& rng, int n) {
  vector<uint16_t> q = independent_queries(rng, n);
  vector<Flat> out;
  for (int leaf = 0; leaf < 8; ++leaf) {
    const int rootbit = (leaf >> 2) & 1;
    const int midbit = (leaf >> 1) & 1;
    const int endbit = leaf & 1;
    int midnode = rootbit ? 2 : 1;
    int endnode;
    if (!rootbit) endnode = midbit ? 4 : 3;
    else endnode = midbit ? 6 : 5;
    out.push_back(make_flat({q[0], q[midnode], q[endnode]},
                            {static_cast<unsigned char>(rootbit), static_cast<unsigned char>(midbit), static_cast<unsigned char>(endbit)}, n));
  }
  return out;
}

static Flat mutate_one_equation(const Flat& old, std::mt19937_64& rng, int n) {
  // Random GL(3,2) row operations first (same flat), then replace one equation.
  array<uint16_t,3> a = old.normal;
  array<unsigned char,3> b = old.rhs;
  for (int s = 0; s < 7; ++s) {
    int i = int(rng() % 3), j = int(rng() % 3); if (i == j) continue;
    a[i] ^= a[j]; b[i] ^= b[j];
  }
  int i = int(rng() % 3);
  const uint16_t lim = uint16_t(1u << n);
  do { a[i] = uint16_t(1 + rng() % (lim - 1)); }
  while (gf2rank({a[0],a[1],a[2]},n) != 3);
  b[i] = static_cast<unsigned char>(rng() & 1);
  return make_flat(a, b, n);
}

static std::string flat_json(const Flat& f, int n) {
  std::ostringstream s;
  s << "{\"normal_rows\":[" << f.normal[0] << ',' << f.normal[1] << ',' << f.normal[2]
    << "],\"rhs\":[" << int(f.rhs[0]) << ',' << int(f.rhs[1]) << ',' << int(f.rhs[2]) << "],\"point_mask_words_hex\":[";
  int words = (1 << n) / 64;
  for (int w = 0; w < words; ++w) { if (w) s << ','; s << "\"0x" << std::hex << std::setw(16) << std::setfill('0') << f.point[w] << std::dec << "\""; }
  s << "]}";
  return s.str();
}

static std::string packing_json(const vector<Flat>& p, int n) {
  std::ostringstream s; s << '[';
  for (int i = 0; i < (int)p.size(); ++i) { if (i) s << ','; s << flat_json(p[i], n); }
  s << ']';
  return s.str();
}

int main(int argc, char** argv) {
  int n = 8, target = 8;
  uint64_t iterations = 1000000, seed = 2026111708ULL;
  std::string output = "packing-pilot.json";
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto take = [&](auto& v) { if (++i >= argc) throw std::runtime_error("missing value"); std::istringstream(argv[i]) >> v; };
    if (a == "--n") take(n); else if (a == "--iterations") take(iterations);
    else if (a == "--seed") take(seed); else if (a == "--target") take(target); else if (a == "--out") take(output);
    else throw std::runtime_error("unknown option " + a);
  }
  if (n < 8 || n > 9 || target < 1 || target > n) throw std::runtime_error("supported: n=8 or 9");
  std::mt19937_64 rng(seed);
  vector<Flat> current = decision_tree_partition(rng, n);
  // Puncture four leaves initially. It is deliberately not the old standard parallel-rank-3 seed.
  std::shuffle(current.begin(), current.end(), rng); current.resize(4);
  vector<Flat> best = current, hit;
  int best_rank = packing_rank(current,n);
  auto start = std::chrono::steady_clock::now();
  uint64_t accepted = 0, restarts = 0, adds = 0, replaces = 0, proposals = 0;
  uint64_t executed = 0;
  for (uint64_t it = 0; it < iterations; ++it) {
    executed = it + 1;
    if (it && it % 50000 == 0) {
      current = decision_tree_partition(rng,n);
      std::shuffle(current.begin(), current.end(), rng);
      current.resize(4 + (rng() % 4)); // 4--7 cells, i.e. puncture 1--4
      ++restarts;
    }
    vector<Flat> cand = current;
    const uint64_t roll = rng() % 100;
    bool moved = false;
    if (roll < 40 && cand.size() < 8) { // insert into presently uncovered points
      Flat f = random_flat(rng,n); ++proposals;
      if (disjoint(f,cand)) { cand.push_back(f); moved = true; ++adds; }
    } else { // replace after deleting one or two cells, allowing rank barriers to be crossed
      int drop = 1 + ((roll < 58 && cand.size() >= 2) ? 1 : 0);
      if ((int)cand.size() > drop) {
        for (int d = 0; d < drop; ++d) cand.erase(cand.begin() + (rng() % cand.size()));
        Flat f = ((rng() & 1) && !current.empty()) ? mutate_one_equation(current[rng()%current.size()],rng,n) : random_flat(rng,n);
        ++proposals;
        if (disjoint(f,cand)) { cand.push_back(f); moved = true; ++replaces; }
      }
    }
    if (!moved) continue;
    const int cr = packing_rank(cand,n), oldr = packing_rank(current,n);
    // Bounded simulated annealing: favour rank, retain a little cardinality and accept finite noise.
    int delta = 100 * (cr - oldr) + int(cand.size()) - int(current.size());
    double temp = 8.0 * (1.0 - double(it) / double(iterations)) + 0.5;
    bool accept = delta >= 0 || std::generate_canonical<double,53>(rng) < std::exp(double(delta) / temp);
    if (accept) { current = std::move(cand); ++accepted; }
    int now = packing_rank(current,n);
    if (now > best_rank || (now == best_rank && current.size() > best.size())) { best = current; best_rank = now; }
    if (now >= target) { hit = current; break; }
  }
  // Exact final recheck: all chosen masks pairwise disjoint and all flats have 2^(n-3) points.
  const vector<Flat>& report = hit.empty() ? best : hit;
  bool valid = true; for (int i=0;i<(int)report.size();++i) {
    int pc=0; for(uint64_t z:report[i].point) pc+=__builtin_popcountll(z);
    if(pc != (1 << (n-3)) || !disjoint(report[i],report,i)) valid=false;
  }
  auto stop = std::chrono::steady_clock::now();
  double seconds = std::chrono::duration<double>(stop-start).count();
  std::ofstream os(output);
  os << "{\n\"status\":\"" << (valid?"passed":"failed") << "\",\n\"method\":\"punctured_depth3_decision_tree_7independent_normals_bounded_anneal\",\n"
     << "\"n\":" << n << ",\"codimension\":3,\"target_normal_rank\":" << target << ",\"seed\":" << seed
     << ",\"iterations_requested\":" << iterations << ",\"iterations_executed\":" << executed
     << ",\"counterexample_found\":" << (!hit.empty()?"true":"false") << ",\"best_normal_rank\":" << packing_rank(report,n)
     << ",\"best_flat_count\":" << report.size() << ",\"pairwise_disjoint_exact\":" << (valid?"true":"false")
     << ",\"accepted_moves\":" << accepted << ",\"add_moves\":" << adds << ",\"replace_moves\":" << replaces
     << ",\"proposals\":" << proposals << ",\"restarts\":" << restarts << ",\"elapsed_seconds\":" << std::setprecision(12) << seconds
     << ",\"flats\":" << packing_json(report,n) << "\n}\n";
  std::cout << "{\"status\":\"" << (valid?"passed":"failed") << "\",\"counterexample_found\":" << (!hit.empty()?"true":"false")
            << ",\"best_normal_rank\":" << packing_rank(report,n) << ",\"best_flat_count\":" << report.size()
            << ",\"elapsed_seconds\":" << seconds << "}\n";
  return valid ? 0 : 2;
}
