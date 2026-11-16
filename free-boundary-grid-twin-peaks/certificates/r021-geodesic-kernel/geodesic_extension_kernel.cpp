// Exact completion-kernel audit for all ordered geometric prefixes on P_n square.
// n is restricted to 3..5: the reverse table has 2^(n^2) unsigned-128 entries.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using U = unsigned __int128;
int n, N;
uint32_t full;
std::vector<uint32_t> adj;
std::vector<U> h;
std::vector<U> numer, denom;
std::vector<uint64_t> prefix_count;
std::vector<int> path;
std::ofstream local_out;

std::string dec(U x) {
  if (!x) return "0";
  std::string s;
  while (x) { s.push_back(char('0' + x % 10)); x /= 10; }
  std::reverse(s.begin(), s.end()); return s;
}
void extend(int root, int last, uint32_t S, int d) {
  prefix_count[d]++;
  denom[d] += h[S];
  uint32_t allowed = adj[last] & ~S;
  U local_numer = 0;
  for (uint32_t b = allowed; b; b &= b - 1) {
    int v = __builtin_ctz(b);
    int dist = std::abs(v % n - root % n) + std::abs(v / n - root / n);
    if (dist != d) continue;
    local_numer += h[S | (1u << v)];
    numer[d] += h[S | (1u << v)];
  }
  local_out << d << '\t';
  for (size_t k=0; k<path.size(); ++k) { if (k) local_out << ','; local_out << path[k]; }
  local_out << '\t' << dec(local_numer) << '\t' << dec(h[S]) << '\n';
  for (uint32_t b = allowed; b; b &= b - 1) {
    int v = __builtin_ctz(b);
    int dist = std::abs(v % n - root % n) + std::abs(v / n - root / n);
    if (dist != d) continue;
    path.push_back(v); extend(root, v, S | (1u << v), d + 1); path.pop_back();
  }
}
int main(int argc, char** argv) {
  if (argc != 3) return 2;
  n = std::atoi(argv[1]); if (n < 3 || n > 5) return 2;
  N = n*n; full = (1u << N) - 1; adj.assign(N, 0);
  for (int v=0; v<N; ++v) for (int u=0;u<N;++u)
    if (std::abs(v%n-u%n)+std::abs(v/n-u/n)==1) adj[v] |= 1u<<u;
  h.assign(size_t(full)+1, 0); h[full]=1;
  for (uint32_t S=full-1;S;--S) {
    uint32_t boundary=0;
    for (uint32_t b=S;b;b&=b-1) boundary |= adj[__builtin_ctz(b)];
    for (uint32_t b=boundary&~S;b;b&=b-1) h[S]+=h[S|(1u<<__builtin_ctz(b))];
  }
  numer.assign(2*n,0); denom.assign(2*n,0); prefix_count.assign(2*n,0);
  local_out.open(argv[2]); if (!local_out) return 3;
  for (int root=0;root<N;++root) {path={root};extend(root,root,1u<<root,1);}
  std::cout << "{\"n\":" << n << ",\"A1\":\"";
  U A1=0;for(int v=0;v<N;v++)A1+=h[1u<<v];std::cout<<dec(A1)<<"\",\"rows\":[";
  for(int d=1;d<2*n-1;d++){if(d>1)std::cout<<',';
    std::cout<<"{\"prefix_vertices\":"<<d<<",\"ordered_prefixes\":"<<prefix_count[d]
      <<",\"prefix_completion_weight\":\""<<dec(denom[d])<<"\",\"extension_completion_weight\":\""<<dec(numer[d])<<"\"}";
  }
  std::cout<<"]}\n";
}
