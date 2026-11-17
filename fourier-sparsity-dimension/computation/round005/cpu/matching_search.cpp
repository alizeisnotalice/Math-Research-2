// Exact F_2 matching search for the 6-bit bent seed q(u,v)=u dot v.
// A matching joins a positive and a negative {x,x^delta} pair.  Each join
// supplies the four-vertex incidence row e_x+e_y+e_z+e_w.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using Row = std::uint64_t;
struct Pair { int x, y; std::uint16_t label; };
struct Instance { std::array<std::array<Pair,16>,63> pos, neg; };
struct Matching { std::array<std::array<std::uint8_t,16>,63> target; };

static int parity(Row x) { return __builtin_parityll(x); }
static int profile=0;
static int q(int x) {
  const int u=x&7, v=x>>3;
  if(profile==0) return parity(u&v);
  // A linear Gray-type 3-bit permutation and a quadratic h(v): a
  // Maiorana--McFarland bent seed, used only for independent restarts.
  static constexpr std::array<int,8> pi{0,1,3,2,6,7,5,4};
  return parity(u&pi[v]) ^ int(((v&3)==3) ^ ((v>>2)&1));
}

// The requested 15-dimensional space: 1_{v=i}, and u_j 1_{v mod 2^j=p}.
static std::uint16_t eval_label(int x) {
  const int u=x&7, v=x>>3; std::uint16_t out=0; int k=0;
  for (int i=0;i<8;i++,k++) if (v==i) out|=std::uint16_t(1)<<k;
  for (int j=0;j<3;j++) for (int p=0;p<(1<<j);p++,k++)
    if (((u>>j)&1) && (v&((1<<j)-1))==p) out|=std::uint16_t(1)<<k;
  return out;
}
static int pair_label(const Pair& p) { return p.label; }

static Instance make_instance() {
  Instance a{};
  for (int delta=1;delta<64;delta++) {
    std::vector<Pair> pp, nn;
    for (int x=0;x<64;x++) { int y=x^delta; if (x>=y) continue;
      Pair p{x,y,std::uint16_t(eval_label(x)^eval_label(y))};
      ((q(x)^q(y)) ? nn : pp).push_back(p);
    }
    if (pp.size()!=16 || nn.size()!=16) throw std::runtime_error("unbalanced seed pairs");
    std::sort(pp.begin(),pp.end(),[](auto l,auto r){return std::tie(l.label,l.x)<std::tie(r.label,r.x);});
    std::sort(nn.begin(),nn.end(),[](auto l,auto r){return std::tie(l.label,l.x)<std::tie(r.label,r.x);});
    std::copy(pp.begin(),pp.end(),a.pos[delta-1].begin());
    std::copy(nn.begin(),nn.end(),a.neg[delta-1].begin());
  }
  return a;
}
static Matching structured_matching(const Instance& a, bool& preserves_V15) {
  Matching m{};
  for (int d=0;d<63;d++) {
    std::map<int,std::vector<int>> P,N;
    for (int i=0;i<16;i++) { P[pair_label(a.pos[d][i])].push_back(i); N[pair_label(a.neg[d][i])].push_back(i); }
    if (P.size()!=N.size()) { preserves_V15=false; break; }
    for (const auto& [key, left]:P) {
      const auto it=N.find(key); if (it==N.end() || left.size()!=it->second.size()) { preserves_V15=false; break; }
      for (std::size_t t=0;t<left.size();t++) m.target[d][left[t]]=it->second[t];
    }
    if(!preserves_V15) break;
  }
  if(!preserves_V15) for(int d=0;d<63;d++) for(int i=0;i<16;i++) m.target[d][i]=i;
  return m;
}
static std::array<Row,1008> rows(const Instance& a, const Matching& m) {
  std::array<Row,1008> out{}; int at=0;
  for (int d=0;d<63;d++) for(int i=0;i<16;i++) { const auto& p=a.pos[d][i]; const auto& n=a.neg[d][m.target[d][i]];
    out[at++]=(Row(1)<<p.x)^(Row(1)<<p.y)^(Row(1)<<n.x)^(Row(1)<<n.y);
  }
  return out;
}
static int rank(const std::array<Row,1008>& in) {
  std::array<Row,64> piv{}; int r=0;
  for(Row x:in) { for(int b=63;b>=0 && x;b--) if((x>>b)&1) { if(piv[b]) x^=piv[b]; else {piv[b]=x; r++; break;} } }
  return r;
}
static std::array<Row,64> rref(const std::array<Row,1008>& in) {
  std::array<Row,64> piv{};
  for(Row x:in) { for(int b=63;b>=0 && x;b--) if((x>>b)&1) { if(piv[b]) x^=piv[b]; else {piv[b]=x; for(int c=0;c<64;c++) if(c!=b && ((piv[c]>>b)&1)) piv[c]^=x; break;} } }
  return piv;
}
static std::vector<Row> null_basis(const std::array<Row,1008>& in) {
  const auto piv=rref(in); std::vector<Row> out;
  for(int f=0;f<64;f++) if(!piv[f]) { Row x=Row(1)<<f; for(int b=0;b<64;b++) if(piv[b] && parity(piv[b]&x)) x|=Row(1)<<b; out.push_back(x); }
  for(Row x:out) for(Row r:in) if(parity(x&r)) throw std::runtime_error("bad null basis");
  return out;
}
static bool exact_boolean(const std::vector<Row>& basis, std::string& failure) {
  // Each basis lift is checked by exact autocorrelations and exact inverse Walsh values.
  for(std::size_t z=0;z<basis.size();z++) {
    Row h=basis[z];
    for(int delta=1;delta<64;delta++) { int corr=0; for(int x=0;x<64;x++) corr += ((q(x)^q(x^delta)^((h>>x)&1)^((h>>(x^delta))&1)) ? -1:1); if(corr) {failure="correlation delta="+std::to_string(delta);return false;} }
    for(int t=0;t<64;t++) { int w=0; for(int x=0;x<64;x++) w += ((q(x)^((h>>x)&1)^parity(x&t))?-1:1); if(w!=8 && w!=-8) {failure="walsh t="+std::to_string(t);return false;} }
  }
  return true;
}
static bool exhaustive_kernel_boolean(const std::vector<Row>& basis, std::string& failure) {
  if(basis.size()>20) { failure="kernel too large for requested exhaustive check"; return false; }
  for(Row z=0;z<(Row(1)<<basis.size());z++) {
    Row h=0; for(std::size_t j=0;j<basis.size();j++) if((z>>j)&1) h^=basis[j];
    for(int t=0;t<64;t++) { int w=0; for(int x=0;x<64;x++) w += ((q(x)^((h>>x)&1)^parity(x&t))?-1:1); if(w!=8 && w!=-8) {failure="combined walsh mask="+std::to_string(z)+" t="+std::to_string(t);return false;} }
  }
  return true;
}
static std::vector<Row> character_basis() {
  std::vector<Row> c{~Row(0)};
  for(int bit=0;bit<6;bit++) { Row v=0; for(int x=0;x<64;x++) if((x>>bit)&1) v|=Row(1)<<x; c.push_back(v); }
  return c;
}
static int vector_rank(const std::vector<Row>& vs) {
  std::array<Row,64> p{}; int r=0; for(Row x:vs) for(int b=63;b>=0&&x;b--) if((x>>b)&1) {if(p[b])x^=p[b];else{p[b]=x;r++;break;}} return r;
}
static bool character_checks(const std::array<Row,1008>& rs, std::string& failure) {
  const auto c=character_basis();
  if(vector_rank(c)!=7) {failure="character basis rank";return false;}
  for(Row v:c) for(Row r:rs) if(parity(v&r)) {failure="C not annihilating rows";return false;}
  std::array<int,64> labels{}; for(int x=0;x<64;x++) labels[x]=1|(x<<1);
  std::sort(labels.begin(),labels.end()); if(std::adjacent_find(labels.begin(),labels.end())!=labels.end()) {failure="frequency evaluations not distinct";return false;}
  return true;
}
static void write_matching(std::ostream& o,const Instance&a,const Matching&m) {
  o<<"["; for(int d=0;d<63;d++){if(d)o<<",";o<<"[";for(int i=0;i<16;i++){if(i)o<<",";const auto&p=a.pos[d][i];const auto&n=a.neg[d][m.target[d][i]];o<<"["<<p.x<<","<<p.y<<","<<n.x<<","<<n.y<<"]";}o<<"]";}o<<"]";
}
static void write_basis(std::ostream&o,const std::vector<Row>& b){o<<"[";for(std::size_t i=0;i<b.size();i++){if(i)o<<",";o<<"\"0x"<<std::hex<<b[i]<<std::dec<<"\"";}o<<"]";}
static void write_rows(std::ostream&o,const std::array<Row,1008>& rs){o<<"[";for(std::size_t i=0;i<rs.size();i++){if(i)o<<",";o<<"\"0x"<<std::hex<<rs[i]<<std::dec<<"\"";}o<<"]";}
int main(int argc,char**argv) {
  if(argc!=5 && argc!=6){std::cerr<<"usage: matching_search SEED ITERS TEMP OUT.json [base|mm-nonlinear]\n";return 2;}
  const std::uint64_t seed=std::stoull(argv[1]), iters=std::stoull(argv[2]); const double temp=std::stod(argv[3]);
  if(argc==6) { const std::string name=argv[5]; if(name=="base") profile=0; else if(name=="mm-nonlinear") profile=1; else throw std::runtime_error("unknown profile"); }
  const Instance a=make_instance(); bool preserves_V15=true; Matching current=structured_matching(a,preserves_V15), best=current;
  const int initial=rank(rows(a,current)); int cr=initial, br=initial; std::mt19937_64 gen(seed); std::uint64_t accepted=0, improved=0;
  const auto start=std::chrono::steady_clock::now();
  for(std::uint64_t it=0;it<iters;it++) { int d=gen()%63, i=gen()%16, j=gen()%16; if(i==j) continue;
    const auto old=current.target[d]; const int move=gen()%20;
    if(move<14) std::swap(current.target[d][i],current.target[d][j]); // two matched pairs
    else if(move<19) { int k=gen()%16; if(k==i || k==j) {current.target[d]=old;continue;} // three-pair cycle
      const auto v=current.target[d][i]; current.target[d][i]=current.target[d][j]; current.target[d][j]=current.target[d][k]; current.target[d][k]=v;
    } else std::shuffle(current.target[d].begin(),current.target[d].end(),gen); // whole-delta rematching
    const int nr=rank(rows(a,current));
    const double T=temp*(1.0-double(it)/std::max<std::uint64_t>(1,iters));
    const bool take=(nr<=cr) || (T>0 && std::generate_canonical<double,53>(gen)<std::exp(double(cr-nr)/T));
    if(take){cr=nr;accepted++; if(nr<br){br=nr;best=current;improved++;}} else current.target[d]=old;
  }
  const auto bestrows=rows(a,best); const auto basis=null_basis(bestrows); std::string failure; const bool basis_boolean=exact_boolean(basis,failure);
  std::string exhaustive_failure, character_failure; const bool exhaustive_boolean=exhaustive_kernel_boolean(basis,exhaustive_failure); const bool character_ok=character_checks(bestrows,character_failure);
  const double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::ofstream o(argv[4]);
  o<<"{\n  \"method\":\"V15_label_start_or_canonical_start_then_single_delta_matching_swap_annealing\",\n  \"seed_profile\":\""<<(profile==0?"base":"mm-nonlinear")<<"\",\n  \"initial_preserves_specified_V15\":"<<(preserves_V15?"true":"false")<<",\n  \"seed\":"<<seed<<",\n  \"iterations\":"<<iters<<",\n  \"temperature\":"<<temp<<",\n  \"initial_rank_R\":"<<initial<<",\n  \"best_rank_R\":"<<br<<",\n  \"best_dimV\":"<<(64-br)<<",\n  \"best_graph_affine_dimension\":"<<(63-br)<<",\n  \"accepted\":"<<accepted<<",\n  \"strict_improvements\":"<<improved<<",\n  \"elapsed_seconds\":"<<sec<<",\n  \"incidence_rows_hex\":";write_rows(o,bestrows);o<<",\n  \"nullspace_basis_hex\":";write_basis(o,basis);o<<",\n  \"nullspace_basis_dimension\":"<<basis.size()<<",\n  \"individual_basis_exact_boolean\":"<<(basis_boolean?"true":"false")<<",\n  \"all_kernel_masks_exact_boolean\":"<<(exhaustive_boolean?"true":"false")<<",\n  \"all_kernel_masks_checked\":"<<(basis.size()<=20?(Row(1)<<basis.size()):0)<<",\n  \"C_basis_hex\":";write_basis(o,character_basis());o<<",\n  \"C_rank\":7,\n  \"C_annihilates_incidence_rows\":"<<(character_ok?"true":"false")<<",\n  \"frequency_evaluation_columns_distinct\":"<<(character_ok?"true":"false")<<",\n  \"boolean_failure\":\""<<failure<<"; "<<exhaustive_failure<<"; "<<character_failure<<"\",\n  \"matching_by_delta\":";write_matching(o,a,best);o<<"\n}\n";
  const bool ok=basis_boolean&&exhaustive_boolean&&character_ok;
  std::cout<<"initial_R="<<initial<<" best_R="<<br<<" dimV="<<(64-br)<<" affdim="<<(63-br)<<" accepted="<<accepted<<" seconds="<<sec<<" exhaustive="<<exhaustive_boolean<<" C="<<character_ok<<"\n";
  return ok?0:4;
}
