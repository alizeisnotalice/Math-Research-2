// Exact subset-DP verifier for finite bipartite prefix identities.
// Remote-only source: this file has deliberately not been compiled or run here.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using U64 = std::uint64_t;
using I64 = std::int64_t;

struct Config {
  std::string mode, out_dir, shard_id, source_sha;
  int left_size = -1;
  U64 start = 0, stop = 0, checkpoint_every = 1000, edge_mask = 0;
};

static std::string now_utc() {
  const auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm tm{};
#if defined(__APPLE__) || defined(__linux__)
  gmtime_r(&t, &tm);
#else
  tm = *std::gmtime(&t);
#endif
  std::ostringstream os; os << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ"); return os.str();
}

static void atomic_write(const fs::path& dest, const std::string& text) {
  fs::create_directories(dest.parent_path());
  // One worker owns each destination checkpoint, so a fixed sibling temp name
  // is sufficient and keeps this source portable without platform PID headers.
  const fs::path tmp = dest.string() + ".tmp";
  { std::ofstream f(tmp); if (!f) throw std::runtime_error("cannot open checkpoint temp"); f << text << '\n'; }
  fs::rename(tmp, dest);
}

static Config parse(int argc, char** argv) {
  Config c;
  std::map<std::string, std::string*> str{{"--mode", &c.mode}, {"--out-dir", &c.out_dir},
      {"--shard-id", &c.shard_id}, {"--source-sha", &c.source_sha}};
  for (int i = 1; i < argc; i += 2) {
    if (i + 1 >= argc) throw std::runtime_error("every option requires a value");
    std::string k(argv[i]), v(argv[i + 1]);
    if (auto it = str.find(k); it != str.end()) *it->second = v;
    else if (k == "--left-size") c.left_size = std::stoi(v);
    else if (k == "--start") c.start = std::stoull(v);
    else if (k == "--stop") c.stop = std::stoull(v);
    else if (k == "--edge-mask") c.edge_mask = std::stoull(v);
    else if (k == "--checkpoint-every") c.checkpoint_every = std::stoull(v);
    else throw std::runtime_error("unknown option: " + k);
  }
  if (c.mode.empty() || c.out_dir.empty() || c.source_sha.empty()) throw std::runtime_error("--mode, --out-dir, --source-sha are required");
  return c;
}

static std::vector<std::vector<int>> masks_by_size(int n) {
  std::vector<std::vector<int>> out(n + 1);
  for (int m = 0; m < (1 << n); ++m) out[__builtin_popcount(static_cast<unsigned>(m))].push_back(m);
  return out;
}

static std::vector<std::pair<int,int>> bipartite_pairs(int n, int left) {
  std::vector<std::pair<int,int>> out;
  for (int x = 0; x < left; ++x) for (int y = left; y < n; ++y) out.emplace_back(x, y);
  return out;
}

static std::vector<int> adjacency(int n, const std::vector<std::pair<int,int>>& pairs, U64 edge_mask) {
  std::vector<int> adj(n, 0);
  for (std::size_t k = 0; k < pairs.size(); ++k) if ((edge_mask >> k) & 1ULL) {
    const auto [x, y] = pairs[k]; adj[x] |= 1 << y; adj[y] |= 1 << x;
  }
  return adj;
}

// h[L] counts valid ordered prefixes rooted at a nonedge (r,p), where the
// first neighbour of p appears at complete rank L+2. h[0] is deliberately
// omitted from the final bins: no first p-neighbour means no short-bin mass.
static std::vector<U64> subset_hist(const std::vector<int>& adj, const std::vector<std::vector<int>>& by_size) {
  const int n = static_cast<int>(adj.size()), all = (1 << n) - 1;
  std::vector<int> neighbor_union(1 << n, 0);
  for (int mask = 1; mask <= all; ++mask) {
    const int bit = mask & -mask;
    neighbor_union[mask] = neighbor_union[mask ^ bit] | adj[__builtin_ctz(static_cast<unsigned>(bit))];
  }
  std::vector<U64> h(n, 0);
  for (int r = 0; r < n; ++r) for (int p = 0; p < n; ++p) {
    if (r == p || (adj[r] & (1 << p))) continue;
    const int initial = (1 << r) | (1 << p);
    std::vector<std::array<U64, 9>> dp(1 << n); // N is at most nine in this program.
    dp[initial][0] = 1;
    for (int rank = 2; rank < n; ++rank) for (int mask : by_size[rank]) {
      const auto row = dp[mask];
      bool nonzero = false; for (int ell = 0; ell <= n - 2; ++ell) nonzero |= row[ell] != 0;
      if (!nonzero) continue;
      int choices = neighbor_union[mask] & ~mask & all;
      while (choices) {
        const int bit = choices & -choices; choices ^= bit;
        const int newmask = mask | bit;
        const bool p_neighbor = (adj[p] & bit) != 0;
        for (int ell = 0; ell <= n - 2; ++ell) if (row[ell]) {
          const int newell = ell ? ell : (p_neighbor ? rank - 1 : 0);
          dp[newmask][newell] += row[ell];
        }
      }
    }
    for (int ell = 1; ell <= n - 2; ++ell) h[ell] += dp[all][ell];
  }
  return h;
}

static std::vector<U64> direct_hist(const std::vector<int>& adj) {
  const int n = static_cast<int>(adj.size());
  std::vector<int> word(n); std::iota(word.begin(), word.end(), 0);
  std::vector<U64> h(n, 0);
  do {
    const int r = word[0], p = word[1];
    if (adj[r] & (1 << p)) continue;
    int seen = (1 << r) | (1 << p), ell = 0; bool valid = true;
    for (int pos = 2; pos < n; ++pos) {
      const int v = word[pos];
      if (!(adj[v] & seen)) { valid = false; break; }
      if (!ell && (adj[p] & (1 << v))) ell = pos - 1;
      seen |= 1 << v;
    }
    if (valid && ell) ++h[ell];
  } while (std::next_permutation(word.begin(), word.end()));
  return h;
}

static I64 r45(const std::vector<U64>& h) { return static_cast<I64>(h[3]) - 3*static_cast<I64>(h[4]) + 2*static_cast<I64>(h[5]); }
static I64 r27(const std::vector<U64>& h) { return static_cast<I64>(h[2]) - 5*static_cast<I64>(h[4]) + 5*static_cast<I64>(h[6]) - 2*static_cast<I64>(h[7]); }
static std::string h_json(const std::vector<U64>& h) {
  std::ostringstream os; os << "["; for (std::size_t i = 0; i < h.size(); ++i) { if (i) os << ","; os << h[i]; } os << "]"; return os.str();
}
static std::string edges_json(const std::vector<std::pair<int,int>>& pairs, U64 mask) {
  std::ostringstream os; os << "["; bool first = true;
  for (std::size_t i = 0; i < pairs.size(); ++i) if ((mask >> i) & 1ULL) { if (!first) os << ","; first=false; os << "[" << pairs[i].first << "," << pairs[i].second << "]"; }
  os << "]"; return os.str();
}

static std::string shard_json(const Config& c, const std::string& status, U64 next, U64 scanned, I64 max_abs, const std::optional<std::string>& counter = std::nullopt, const std::string& failure = "") {
  std::ostringstream os;
  os << "{\"status\":\"" << status << "\",\"mode\":\"shard-r27\",\"N\":9,\"left_size\":" << c.left_size
     << ",\"range_start\":" << c.start << ",\"range_stop_exclusive\":" << c.stop << ",\"next_mask\":" << next
     << ",\"graphs_scanned\":" << scanned << ",\"max_abs_R27\":" << max_abs << ",\"shard_id\":\"" << c.shard_id
     << "\",\"source_sha256\":\"" << c.source_sha << "\",\"utc\":\"" << now_utc() << "\"";
  if (counter) os << ",\"counterexample\":" << *counter;
  if (!failure.empty()) os << ",\"failure\":\"" << failure << "\"";
  return os << "}", os.str();
}

static void scan_shard(const Config& c) {
  if (c.left_size < 1 || c.left_size > 4 || c.shard_id.empty()) throw std::runtime_error("shard requires left size 1..4 and shard id");
  const auto pairs = bipartite_pairs(9, c.left_size);
  const U64 total = 1ULL << pairs.size();
  if (c.stop == 0 || c.checkpoint_every == 0) throw std::runtime_error("shard requires nonzero --stop and --checkpoint-every");
  if (c.start > c.stop || c.stop > total) throw std::runtime_error("invalid shard range");
  const auto checkpoint = fs::path(c.out_dir) / ("shard-l" + std::to_string(c.left_size) + "-" + c.shard_id + ".json");
  // A resumed suffix must use a new shard id. Never erase a running, failed,
  // completed, or counterexample receipt merely because a launcher restarted.
  if (fs::exists(checkpoint)) throw std::runtime_error("refusing to overwrite existing checkpoint: " + checkpoint.string());
  const auto by_size = masks_by_size(9); U64 scanned = 0; I64 max_abs = 0;
  try {
    atomic_write(checkpoint, shard_json(c, "running", c.start, 0, 0));
    for (U64 mask = c.start; mask < c.stop; ++mask) {
      const auto h = subset_hist(adjacency(9, pairs, mask), by_size);
      const I64 rel = r27(h); ++scanned; max_abs = std::max(max_abs, rel < 0 ? -rel : rel);
      if (rel != 0) {
        std::ostringstream ce; ce << "{\"edge_mask\":" << mask << ",\"edges\":" << edges_json(pairs, mask)
                                  << ",\"h\":" << h_json(h) << ",\"R27\":" << rel << "}";
        atomic_write(checkpoint, shard_json(c, "counterexample", mask + 1, scanned, max_abs, ce.str()));
        std::cerr << "COUNTEREXAMPLE saved in " << checkpoint << "\n";
        std::exit(3); // Do not proceed past a candidate counterexample.
      }
      if (scanned % c.checkpoint_every == 0) atomic_write(checkpoint, shard_json(c, "running", mask + 1, scanned, max_abs));
    }
    atomic_write(checkpoint, shard_json(c, "completed", c.stop, scanned, max_abs));
  } catch (const std::exception& e) {
    atomic_write(checkpoint, shard_json(c, "failed", c.start + scanned, scanned, max_abs, std::nullopt, e.what()));
    throw;
  }
}

static std::vector<std::pair<int,int>> sample_edges(const std::string& name) {
  std::vector<std::pair<int,int>> e;
  if (name == "P9") for (int i=0;i<8;++i) e.emplace_back(i,i+1);
  else if (name == "K23-tail") { for (int x=0;x<2;++x) for (int y=2;y<5;++y) e.emplace_back(x,y); e.insert(e.end(), {{2,5},{5,6},{6,7},{7,8}}); }
  else if (name == "C6-tail") { for (int i=0;i<6;++i) e.emplace_back(i,(i+1)%6); e.insert(e.end(), {{0,6},{6,7},{7,8}}); }
  else throw std::runtime_error("unknown sample");
  return e;
}

static std::vector<int> sample_adj(const std::string& name) {
  const auto e = sample_edges(name); std::vector<int> a(9, 0);
  for (const auto [x,y] : e) { a[x] |= 1<<y; a[y] |= 1<<x; }
  return a;
}

static void gate_n7(const Config& c) {
  std::ostringstream out; out << "{\"status\":\"passed\",\"mode\":\"gate-n7-r45\",\"source_sha256\":\"" << c.source_sha << "\",\"families\":[";
  bool first = true; U64 scanned = 0;
  for (int left=1; left<=3; ++left) {
    const auto pairs = bipartite_pairs(7,left);
    const auto by_size = masks_by_size(7);
    const U64 total=1ULL<<pairs.size();
    for (U64 mask=0;mask<total;++mask) {
      const auto h = subset_hist(adjacency(7,pairs,mask),by_size); const I64 rel = r45(h);
      if (rel != 0) {
        std::ostringstream ce; ce << "{\"status\":\"counterexample\",\"mode\":\"gate-n7-r45\",\"left_size\":"<<left<<",\"edge_mask\":"<<mask<<",\"edges\":"<<edges_json(pairs,mask)<<",\"h\":"<<h_json(h)<<",\"R45\":"<<rel<<",\"source_sha256\":\""<<c.source_sha<<"\",\"utc\":\""<<now_utc()<<"\"}";
        atomic_write(fs::path(c.out_dir)/"gate-n7-r45-counterexample.json",ce.str()); std::exit(3);
      }
    }
    if (!first) out << ","; first=false; out << "{\"left_size\":"<<left<<",\"edge_subsets\":"<<total<<",\"nonzero_R45\":0}"; scanned += total;
  }
  out << "],\"graphs_scanned\":" << scanned << ",\"utc\":\""<<now_utc()<<"\"}";
  atomic_write(fs::path(c.out_dir)/"gate-n7-r45.json",out.str());
}

static void gate_n9(const Config& c) {
  std::ostringstream out; out << "{\"status\":\"passed\",\"mode\":\"gate-n9-direct\",\"source_sha256\":\""<<c.source_sha<<"\",\"samples\":[";
  bool first=true;
  for (const std::string name : {"P9","K23-tail","C6-tail"}) {
    const auto e = sample_edges(name);
    const auto a = sample_adj(name);
    const auto dp=subset_hist(a,masks_by_size(9)); const auto direct=direct_hist(a);
    if (dp != direct) throw std::runtime_error("subset DP/direct 9! mismatch for "+name);
    const I64 rel=r27(dp);
    if (!first) out << ","; first=false;
    out << "{\"name\":\""<<name<<"\",\"edges\":"<<edges_json(e,(1ULL<<e.size())-1)<<",\"h\":"<<h_json(dp)<<",\"R27\":"<<rel<<",\"dp_equals_direct\":true}";
    if (rel != 0) { out << "],\"counterexample_sample\":\""<<name<<"\",\"utc\":\""<<now_utc()<<"\"}"; atomic_write(fs::path(c.out_dir)/"gate-n9-direct-counterexample.json",out.str()); std::exit(3); }
  }
  out << "],\"utc\":\""<<now_utc()<<"\"}"; atomic_write(fs::path(c.out_dir)/"gate-n9-direct.json",out.str());
}

static void verify_edge_mask(const Config& c) {
  if (c.left_size < 1 || c.left_size > 4) throw std::runtime_error("verify mode requires --left-size 1..4");
  const auto pairs = bipartite_pairs(9, c.left_size);
  if (c.edge_mask >= (1ULL << pairs.size())) throw std::runtime_error("edge mask is outside host range");
  const auto adj = adjacency(9, pairs, c.edge_mask);
  const auto dp = subset_hist(adj, masks_by_size(9));
  const auto direct = direct_hist(adj);
  if (dp != direct) throw std::runtime_error("independent direct 9! check disagrees with subset DP");
  const I64 rel = r27(dp);
  std::ostringstream out;
  out << "{\"status\":\"" << (rel ? "confirmed-counterexample" : "verified-zero")
      << "\",\"mode\":\"verify-edge-mask-r27\",\"left_size\":" << c.left_size
      << ",\"edge_mask\":" << c.edge_mask << ",\"edges\":" << edges_json(pairs,c.edge_mask)
      << ",\"h\":" << h_json(dp) << ",\"R27\":" << rel << ",\"dp_equals_direct\":true"
      << ",\"source_sha256\":\"" << c.source_sha << "\",\"utc\":\"" << now_utc() << "\"}";
  atomic_write(fs::path(c.out_dir) / ("direct-verify-l" + std::to_string(c.left_size) + "-mask" + std::to_string(c.edge_mask) + ".json"), out.str());
}

int main(int argc, char** argv) {
  try { const Config c=parse(argc,argv); if(c.mode=="shard-r27") scan_shard(c); else if(c.mode=="gate-n7-r45") gate_n7(c); else if(c.mode=="gate-n9-direct") gate_n9(c); else if(c.mode=="verify-edge-mask-r27") verify_edge_mask(c); else throw std::runtime_error("mode must be shard-r27, gate-n7-r45, gate-n9-direct, or verify-edge-mask-r27"); return 0; }
  catch(const std::exception& e) { std::cerr << "fatal: " << e.what() << "\n"; return 2; }
}
