// Low-memory one-peak Psi search with bounded local proposal rollback.
// Search arithmetic is double; verify_psi.py independently checks the final word.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct State {
  int n, N, root;
  std::vector<int> a, pos, l, b, sole;
  double psi = 0;
  std::vector<std::vector<int>> nb;

  struct Undo {
    int i = -1, j = -1;
    double psi = 0;
    std::vector<int> S, P, sole0, l0, b0;
  };

  State(int q, int r)
      : n(q), N(q * q), root(r), a(N), pos(N), l(N), b(N), sole(N), nb(N) {
    for (int v = 0; v < N; ++v) {
      int x = v % n, y = v / n;
      for (auto d : {std::pair{-1, 0}, std::pair{1, 0}, std::pair{0, -1},
                     std::pair{0, 1}}) {
        int X = x + d.first, Y = y + d.second;
        if (X >= 0 && X < n && Y >= 0 && Y < n) nb[v].push_back(Y * n + X);
      }
    }
  }

  bool eligible(int p) const {
    return p != root && std::find(nb[root].begin(), nb[root].end(), p) == nb[root].end();
  }
  int earlycount(int v) const {
    int c = 0;
    for (int u : nb[v]) c += pos[u] < pos[v];
    return c;
  }
  int soleparent(int v) const {
    int z = -1;
    for (int u : nb[v]) {
      if (pos[u] < pos[v]) {
        if (z >= 0) return -1;
        z = u;
      }
    }
    return z;
  }
  void recompute_p(int p) {
    l[p] = N + 1;
    for (int u : nb[p]) l[p] = std::min(l[p], pos[u]);
    b[p] = N + 1;
    for (int v : nb[p]) if (sole[v] == p) b[p] = std::min(b[p], pos[v]);
  }
  double term(int p) const {
    if (!eligible(p)) return 0;
    int R = b[p] - l[p] - 1;
    if (R <= 0) throw std::runtime_error("bad R");
    return double(l[p] - 1) / R;
  }
  void rebuild() {
    for (int v = 0; v < N; ++v) {
      if (v != root && earlycount(v) == 0) throw std::runtime_error("not onepeak");
      sole[v] = soleparent(v);
    }
    psi = 0;
    for (int p = 0; p < N; ++p) {
      recompute_p(p);
      psi += term(p);
    }
  }
  void init_bfs() {
    std::vector<int> d(N, -1);
    d[root] = 0;
    a[0] = root;
    for (int i = 0, k = 1; i < k; ++i) {
      for (int u : nb[a[i]]) if (d[u] < 0) {
        d[u] = d[a[i]] + 1;
        a[k++] = u;
      }
    }
    for (int i = 0; i < N; ++i) pos[a[i]] = i + 1;
    rebuild();
  }
  void init_row() {
    if (root != 0) throw std::runtime_error("row root");
    for (int i = 0; i < N; ++i) a[i] = i, pos[i] = i + 1;
    rebuild();
  }
  void init_ring() {
    if (root != 0) throw std::runtime_error("ring root");
    int z = 0;
    a[z++] = 0;
    for (int k = 1; k < n; ++k) {
      for (int y = 0; y <= k; ++y) a[z++] = y * n + k;
      for (int x = k - 1; x >= 0; --x) a[z++] = k * n + x;
    }
    for (int i = 0; i < N; ++i) pos[a[i]] = i + 1;
    rebuild();
  }

  void init_word(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("input word open");
    std::string header, tag, nr, rr;
    std::getline(f, header);
    std::istringstream hs(header);
    int q = -1, r = -1;
    if (!(hs >> nr >> q >> rr >> r) || nr != "n" || rr != "root" || q != n || r != root)
      throw std::runtime_error("input header mismatch");
    if (!(f >> tag) || tag != "word") throw std::runtime_error("input word tag");
    std::vector<bool> seen(N);
    for (int i = 0; i < N; ++i) {
      int v = -1;
      if (!(f >> v) || v < 0 || v >= N || seen[v]) throw std::runtime_error("input not permutation");
      a[i] = v; pos[v] = i + 1; seen[v] = true;
    }
    if (f >> tag) throw std::runtime_error("input trailing token");
    if (a[0] != root) throw std::runtime_error("input root");
    rebuild();
  }

  // Snapshot exactly the arrays whose values can change under a transposition.
  // S: endpoints and their graph neighbors (sole-parent values); P: S, former
  // sole parents, and one graph-neighborhood expansion (l,b values).
  bool swap_try(int i, int j, Undo& u) {
    if (i == j || i == 0 || j == 0) return false;
    u = Undo{};
    u.i = i;
    u.j = j;
    u.psi = psi;
    u.S = {a[i], a[j]};
    const size_t endpoints = u.S.size();
    for (size_t z = 0; z < endpoints; ++z) {
      for (int x : nb[u.S[z]]) u.S.push_back(x);
    }
    std::sort(u.S.begin(), u.S.end());
    u.S.erase(std::unique(u.S.begin(), u.S.end()), u.S.end());
    u.P = u.S;
    for (int v : u.S) if (sole[v] >= 0) u.P.push_back(sole[v]);
    const size_t parents_before_expansion = u.P.size();
    for (size_t z = 0; z < parents_before_expansion; ++z) {
      for (int x : nb[u.P[z]]) u.P.push_back(x);
    }
    std::sort(u.P.begin(), u.P.end());
    u.P.erase(std::unique(u.P.begin(), u.P.end()), u.P.end());
    u.sole0.reserve(u.S.size());
    for (int v : u.S) u.sole0.push_back(sole[v]);
    u.l0.reserve(u.P.size());
    u.b0.reserve(u.P.size());
    for (int p : u.P) {
      u.l0.push_back(l[p]);
      u.b0.push_back(b[p]);
    }

    double old = 0;
    for (int p : u.P) old += term(p);
    std::swap(a[i], a[j]);
    pos[a[i]] = i + 1;
    pos[a[j]] = j + 1;
    for (int v : u.S) sole[v] = soleparent(v);
    for (int v : u.S) if (v != root && earlycount(v) == 0) return false;
    for (int p : u.P) recompute_p(p);
    double now = 0;
    for (int p : u.P) now += term(p);
    psi += now - old;
    return true;
  }
  void rollback(const Undo& u) {
    if (u.i < 0) return;
    std::swap(a[u.i], a[u.j]);
    pos[a[u.i]] = u.i + 1;
    pos[a[u.j]] = u.j + 1;
    psi = u.psi;
    for (size_t z = 0; z < u.S.size(); ++z) sole[u.S[z]] = u.sole0[z];
    for (size_t z = 0; z < u.P.size(); ++z) {
      l[u.P[z]] = u.l0[z];
      b[u.P[z]] = u.b0[z];
    }
  }
};

int main(int c, char** v) {
  try {
    if (c < 6 || c > 8) throw std::runtime_error("usage n root seed proposals out [audit] [bfs|row|ring]");
    int n = std::stoi(v[1]), r = std::stoi(v[2]);
    uint64_t seed = std::stoull(v[3]), steps = std::stoull(v[4]);
    if (n < 2 || r < 0 || r >= n * n || steps == 0) throw std::runtime_error("bounds");
    int chk = c >= 7 ? std::stoi(v[6]) : 1000;
    if (chk < 1) throw std::runtime_error("audit");
    State s(n, r);
    std::string init = c == 8 ? v[7] : "bfs";
    if (init == "bfs") s.init_bfs();
    else if (init == "row") s.init_row();
    else if (init == "ring") s.init_ring();
    else if (init.rfind("word:", 0) == 0) s.init_word(init.substr(5));
    else throw std::runtime_error("init");
    auto audit = [](State& x) {
      double old = x.psi;
      State z = x;
      z.rebuild();
      double tol = 1e-11 * std::max(1.0, std::abs(z.psi));
      if (std::abs(old - z.psi) > tol) throw std::runtime_error("incremental mismatch");
      if (x.a != z.a || x.pos != z.pos || x.l != z.l || x.b != z.b || x.sole != z.sole) {
        throw std::runtime_error("incremental array mismatch");
      }
      x = std::move(z);
    };
    State best = s;
    std::mt19937_64 g(seed);
    uint64_t acc = 0;
    for (uint64_t z = 0; z < steps; ++z) {
      int i = 1 + g() % (s.N - 1), j = 1 + g() % (s.N - 1);
      double before = s.psi;
      State::Undo undo;
      if (s.swap_try(i, j, undo)) {
        if (z % chk == 0) audit(s);
        double after = s.psi;
        if (after <= before || std::log(std::generate_canonical<double, 53>(g)) < -(after - before) / .05) {
          ++acc;
        } else {
          s.rollback(undo);
        }
      } else {
        s.rollback(undo);
      }
      if (s.psi < best.psi) best = s;
    }
    audit(s);
    audit(best);
    std::ofstream f(v[5]);
    if (!f) throw std::runtime_error("out");
    f.precision(17);
    f << "n " << n << " root " << r << " seed " << seed << " proposals " << steps
      << " accepted " << acc << " psi " << best.psi << " norm "
      << best.psi / (best.N * std::log(double(n))) << " init " << init << "\nword";
    for (int x : best.a) f << ' ' << x;
    f << '\n';
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
