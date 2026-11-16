// Exact enumeration of g(u,y)=sigma(u)(-1)^{b(u).y}, u in F_2^3, y in F_2^2.
// A global sign leaves support unchanged; normalize sigma(0)=+1.
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

static int parity(unsigned x) { return __builtin_parity(x); }
static int f2rank(unsigned mask) {
  unsigned piv[5]{};
  for (unsigned v = 0; v < 32; ++v) if (mask & (1u << v)) {
    unsigned x = v;
    for (int bit = 4; bit >= 0; --bit) if ((x >> bit) & 1u) {
      if (piv[bit]) x ^= piv[bit]; else { piv[bit] = x; break; }
    }
  }
  int rank = 0; for (unsigned x : piv) rank += x != 0; return rank;
}
static std::string choices(uint64_t id) {
  constexpr uint64_t tail_count = 2097152; // 8^7
  std::string out = std::to_string(id / tail_count); id %= tail_count;
  for (int u = 1; u < 8; ++u) { out += "," + std::to_string(id % 8); id /= 8; }
  return out;
}
int main(int argc, char** argv) {
  if (argc != 4) { std::cerr << "usage: selector_search START COUNT OUT.json\n"; return 2; }
  const uint64_t start = std::stoull(argv[1]), count = std::stoull(argv[2]);
  constexpr uint64_t total = 8388608; // 4*8^7
  if (start > total || count > total - start) return 3;
  std::array<uint64_t, 6> seen{}, bestid{};
  std::array<int, 6> bests; bests.fill(std::numeric_limits<int>::max());
  std::array<unsigned, 6> bestmask{};
  const auto begun = std::chrono::steady_clock::now();
  for (uint64_t id = start; id < start + count; ++id) {
    std::array<unsigned, 8> b{}; std::array<int, 8> sigma{};
    uint64_t q = id; b[0] = q / 2097152; sigma[0] = 1; q %= 2097152;
    for (unsigned u = 1; u < 8; ++u) { const unsigned z = q % 8; q /= 8; b[u] = z & 3u; sigma[u] = z < 4 ? 1 : -1; }
    unsigned mask = 0;
    for (unsigned c = 0; c < 4; ++c) for (unsigned a = 0; a < 8; ++a) {
      int walsh = 0;
      for (unsigned u = 0; u < 8; ++u) if (b[u] == c) walsh += sigma[u] * (parity(a & u) ? -1 : 1);
      if (walsh) mask |= 1u << ((a << 2) | c);
    }
    const int sparsity = __builtin_popcount(mask), dimension = f2rank(mask);
    if (dimension) {
      ++seen[dimension];
      if (sparsity < bests[dimension] || (sparsity == bests[dimension] && id < bestid[dimension])) {
        bests[dimension] = sparsity; bestid[dimension] = id; bestmask[dimension] = mask;
      }
    }
  }
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-begun).count();
  std::ofstream out(argv[3]);
  out << "{\n  \"family\": \"signed-selector-u3-y2\",\n  \"global_sign_normalized\": true,\n  \"start\": " << start << ",\n  \"count\": " << count << ",\n  \"total_normalized\": " << total << ",\n  \"elapsed_seconds\": " << seconds << ",\n  \"by_rank\": [";
  for (int d = 1; d <= 5; ++d) {
    if (d > 1) out << ',';
    out << "{\"d\": " << d << ", \"seen\": " << seen[d];
    if (seen[d]) out << ", \"min_s\": " << bests[d] << ", \"id\": " << bestid[d] << ", \"choices\": \"" << choices(bestid[d]) << "\", \"support_mask_hex\": \"0x" << std::hex << bestmask[d] << std::dec << "\"";
    out << '}';
  }
  out << "]\n}\n";
}
