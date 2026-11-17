// Exact F2/ integer-FWHT audit for the R18 non-affine four-block family.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

using U = std::uint64_t;
struct Bits { std::array<U, 4> w{}; };

static bool getbit(const Bits& x, int p) { return (x.w[p >> 6] >> (p & 63)) & 1ULL; }
static void setbit(Bits& x, int p) { x.w[p >> 6] |= 1ULL << (p & 63); }
static void togglebit(Bits& x, int p) { x.w[p >> 6] ^= 1ULL << (p & 63); }
static Bits xorb(Bits x, const Bits& y) { for (int i = 0; i < 4; ++i) x.w[i] ^= y.w[i]; return x; }
static bool zero(const Bits& x) { return x.w == std::array<U,4>{}; }
static int pop(int x) { int n=0; while (x) { n += x&1; x >>= 1; } return n; }

static void fwht(std::vector<int>& a) {
  for (int h = 1; h < static_cast<int>(a.size()); h <<= 1)
    for (int p = 0; p < static_cast<int>(a.size()); p += 2*h)
      for (int i = p; i < p+h; ++i) { int u=a[i], v=a[i+h]; a[i]=u+v; a[i+h]=u-v; }
}

template <class Fn> static void parallel_for(int count, int threads, Fn fn) {
  std::atomic<int> next{0}; std::vector<std::thread> pool;
  for (int t=0;t<threads;++t) pool.emplace_back([&] {
    for (;;) { int start=next.fetch_add(32); if (start>=count) return;
      for (int i=start;i<std::min(count,start+32);++i) fn(i); }
  });
  for (auto& t:pool) t.join();
}

static int f_value(int x, int r) {
  int x0=(x>>0)&1, x1=(x>>1)&1, x2=(x>>2)&1, x3=(x>>3)&1, x4=(x>>4)&1, x5=(x>>5)&1;
  int q=(x1&x2) ^ (x3&x4);
  int f=(x1&x2) ^ (x0&(x1^x3)) ^ ((x5^q)&(x2^x4));
  for (int k=0;k<r;++k) f ^= ((x>>(6+k))&1) & ((x>>(6+r+k))&1);
  return f;
}

static int q_value(int x) { return (((x>>1)&1)&((x>>2)&1)) ^ (((x>>3)&1)&((x>>4)&1)); }

int main(int argc, char** argv) {
  if (argc != 4) { std::cerr << "usage: verify_nonaffine_blocks r(0|1) threads output.json\n"; return 2; }
  const int r=std::atoi(argv[1]), threads=std::atoi(argv[2]);
  if ((r!=0 && r!=1) || threads<1 || threads>15) return 2;
  const int n=6+2*r, N=1<<n, limbs=N/64, expectedD=(r==0?9:14), expectedAbs=1<<(3+r);
  const int blockCount=4*(1<<r), blockSize=16*(1<<r);
  std::vector<Bits> blocks(blockCount), generators;
  Bits F;
  for (int x=0;x<N;++x) {
    if (f_value(x,r)) setbit(F,x);
    int a=(x>>0)&1, b=((x>>5)&1)^q_value(x), c=(r?((x>>(6+r))&1):0);
    int index=(a<<1)|b; index |= c<<2;
    setbit(blocks[index],x);
  }
  // Independent insertion by high pivot bits gives a basis of the requested V.
  std::array<Bits,256> pivot{}; std::array<bool,256> has{};
  auto add_generator=[&](Bits x) {
    for (int p=N-1;p>=0;--p) if (getbit(x,p)) {
      if (has[p]) x=xorb(x,pivot[p]); else { pivot[p]=x; has[p]=true; generators.push_back(x); return; }
    }
  };
  for (const auto& b:blocks) add_generator(b);
  Bits constant; for (int x=0;x<N;++x) setbit(constant,x); add_generator(constant);
  for (int k=0;k<n;++k) { Bits affine; for (int x=0;x<N;++x) if ((x>>k)&1) setbit(affine,x); add_generator(affine); }
  const int D=static_cast<int>(generators.size());
  std::vector<Bits> codewords(1<<D);
  for (int k=0;k<D;++k) for (int s=0;s<(1<<k);++s) codewords[s|(1<<k)]=xorb(codewords[s],generators[k]);
  std::atomic<int> badPhases{0}, badTerms{0}, maximumDeviation{0};
  auto began=std::chrono::steady_clock::now();
  parallel_for(static_cast<int>(codewords.size()),threads,[&](int phase) {
    Bits mask=xorb(F,codewords[phase]); std::vector<int> a(N);
    for (int x=0;x<N;++x) a[x]=getbit(mask,x)?-1:1;
    fwht(a); int bad=0, dev=0;
    for (int v:a) { int d=std::abs(std::abs(v)-expectedAbs); if (d) ++bad; dev=std::max(dev,d); }
    if (bad) { badPhases.fetch_add(1); badTerms.fetch_add(bad); int old=maximumDeviation.load(); while(dev>old && !maximumDeviation.compare_exchange_weak(old,dev)){} }
  });
  // Check dual partial supports.  For r=0 this is direct for every V phase;
  // for r=1 direct for every global linear representative.  Adding one block
  // indicator negates exactly that block's restricted signal, hence cannot
  // change its support; adding a constant also cannot change a support.
  const int partialCases=(r==0?static_cast<int>(codewords.size()):N);
  std::atomic<int> partialFailures{0}, partialBlockChecks{0}, partialSizeFailures{0}, partialMagnitudeFailures{0};
  parallel_for(partialCases,threads,[&](int phase) {
    Bits mask=F;
    if (r==0) mask=xorb(mask,codewords[phase]);
    else for (int x=0;x<N;++x) if (pop(x&phase)&1) togglebit(mask,x);
    std::vector<unsigned char> seen(N,0); bool fail=false;
    for (int k=0;k<blockCount;++k) {
      std::vector<int> a(N,0);
      for (int x=0;x<N;++x) if (getbit(blocks[k],x)) a[x]=getbit(mask,x)?-1:1;
      fwht(a);
      int supportSize=0; bool magnitudeFailure=false;
      for (int xi=0;xi<N;++xi) if (a[xi]) {
        ++supportSize; if (std::abs(a[xi])!=expectedAbs) magnitudeFailure=true;
        if (seen[xi]) fail=true; seen[xi]=1;
      }
      partialBlockChecks.fetch_add(1);
      if (supportSize!=blockSize) { partialSizeFailures.fetch_add(1); fail=true; }
      if (magnitudeFailure) { partialMagnitudeFailures.fetch_add(1); fail=true; }
    }
    for (int xi=0;xi<N;++xi) if (!seen[xi]) fail=true;
    if (fail) partialFailures.fetch_add(1);
  });
  int partitionFailures=0, blockWeightFailures=0;
  for (int x=0;x<N;++x) { int count=0; for (const auto& b:blocks) count+=getbit(b,x); if (count!=1) ++partitionFailures; }
  for (const auto& b:blocks) { int weight=0; for(int x=0;x<N;++x) weight+=getbit(b,x); if(weight!=blockSize) ++blockWeightFailures; }
  bool generator_rank_ok=(D==expectedD), bent_ok=(badPhases.load()==0), partial_ok=(partialFailures.load()==0);
  double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
  std::ofstream out(argv[3]); if(!out) return 3;
  out << "{\n"
      << "  \"schema\": \"a28-r018-nonaffine-block-cpu-v1\",\n"
      << "  \"r\": "<<r<<",\n  \"n\": "<<n<<",\n  \"truth_table_size\": "<<N<<",\n"
      << "  \"block_count\": "<<blockCount<<",\n  \"block_size\": "<<blockSize<<",\n"
      << "  \"generator_rank\": "<<D<<",\n  \"expected_generator_rank\": "<<expectedD<<",\n"
      << "  \"generator_rank_ok\": "<<(generator_rank_ok?"true":"false")<<",\n"
      << "  \"phase_count_expected\": "<<(1<<D)<<",\n  \"phase_count_processed\": "<<codewords.size()<<",\n"
      << "  \"expected_walsh_magnitude\": "<<expectedAbs<<",\n  \"bad_phase_count\": "<<badPhases.load()<<",\n"
      << "  \"bad_walsh_term_count\": "<<badTerms.load()<<",\n  \"maximum_abs_deviation\": "<<maximumDeviation.load()<<",\n"
      << "  \"partial_support_direct_cases\": "<<partialCases<<",\n  \"partial_support_failure_count\": "<<partialFailures.load()<<",\n"
      << "  \"partial_block_checks\": "<<partialBlockChecks.load()<<",\n"
      << "  \"partial_support_size_failure_count\": "<<partialSizeFailures.load()<<",\n"
      << "  \"partial_nonzero_magnitude_failure_count\": "<<partialMagnitudeFailures.load()<<",\n"
      << "  \"block_partition_failure_count\": "<<partitionFailures<<",\n  \"block_weight_failure_count\": "<<blockWeightFailures<<",\n"
      << "  \"elapsed_seconds\": "<<elapsed<<",\n"
      << "  \"all_checks_passed\": "<<(generator_rank_ok&&bent_ok&&partial_ok&&partitionFailures==0&&blockWeightFailures==0?"true":"false")<<"\n}\n";
  out.close();
  return (generator_rank_ok&&bent_ok&&partial_ok&&partitionFailures==0&&blockWeightFailures==0)?0:10;
}
