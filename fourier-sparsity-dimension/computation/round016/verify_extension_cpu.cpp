// Exact 15-thread CPU verifier for the R16 two-block extension family.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using U = std::uint64_t;

static constexpr U TOGGLE_Z = 0xffff0000ffff0000ULL;

static int bit(U mask, int position) { return static_cast<int>((mask >> position) & 1ULL); }

static U f0_six_mask() {
  U out = 0;
  for (int x = 0; x < 64; ++x) {
    int f = (((x >> 0) & 1) & ((x >> 2) & 1)) ^
            (((x >> 1) & 1) & ((x >> 3) & 1)) ^
            (((x >> 4) & 1) & ((x >> 5) & 1));
    out |= static_cast<U>(f) << x;
  }
  return out;
}

static void fwht(std::array<int, 64>& values) {
  for (int h = 1; h < 64; h <<= 1) {
    for (int base = 0; base < 64; base += 2 * h) {
      for (int j = base; j < base + h; ++j) {
        int left = values[j], right = values[j + h];
        values[j] = left + right;
        values[j + h] = left - right;
      }
    }
  }
}

static U min_atomic(std::atomic<U>& target, U value) {
  U old = target.load(std::memory_order_relaxed);
  while (value < old && !target.compare_exchange_weak(old, value, std::memory_order_relaxed)) {}
  return old;
}

static U extension_mask(U h0, U h1, int a) {
  U g = h0 | (h0 << 16) | (h1 << 32) | (h1 << 48);
  return a ? (g ^ TOGGLE_Z) : g;
}

static int check_phase(U function_mask, int* first_bad_frequency, int* first_bad_value,
                       int* maximum_deviation) {
  std::array<int, 64> values{};
  for (int x = 0; x < 64; ++x) values[x] = bit(function_mask, x) ? -1 : 1;
  fwht(values);
  int bad = 0;
  for (int frequency = 0; frequency < 64; ++frequency) {
    int magnitude = std::abs(values[frequency]);
    *maximum_deviation = std::max(*maximum_deviation, std::abs(magnitude - 8));
    if (magnitude != 8) {
      ++bad;
      if (*first_bad_frequency < 0) {
        *first_bad_frequency = frequency;
        *first_bad_value = values[frequency];
      }
    }
  }
  return bad;
}

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "usage: verify_extension_cpu seed_masks.txt threads output.json\n";
    return 2;
  }
  const int threads = std::atoi(argv[2]);
  if (threads < 1 || threads > 15) return 2;
  std::ifstream input(argv[1]);
  std::vector<U> seeds;
  U value;
  while (input >> value) {
    if (value >= (1ULL << 16)) return 3;
    seeds.push_back(value);
  }
  if (!input.eof() || seeds.size() != 896) return 3;
  const U n = seeds.size();
  const U total = 2 * n * n;
  const U base = f0_six_mask();
  std::atomic<U> next{0}, processed{0}, bad_phase_count{0}, bad_term_count{0};
  std::atomic<U> first_bad{std::numeric_limits<U>::max()};
  std::atomic<int> maximum_bad_deviation{0};
  const auto began = std::chrono::steady_clock::now();
  auto worker = [&]() {
    constexpr U chunk = 128;
    while (true) {
      U start = next.fetch_add(chunk, std::memory_order_relaxed);
      if (start >= total) return;
      U stop = std::min(total, start + chunk);
      for (U index = start; index < stop; ++index) {
        U pair = index % (n * n);
        int a = static_cast<int>(index / (n * n));
        U i = pair / n, j = pair % n;
        U function_mask = base ^ extension_mask(seeds[i], seeds[j], a);
        int freq = -1, walsh = 0, deviation = 0;
        int bad = check_phase(function_mask, &freq, &walsh, &deviation);
        if (bad) {
          bad_phase_count.fetch_add(1, std::memory_order_relaxed);
          bad_term_count.fetch_add(static_cast<U>(bad), std::memory_order_relaxed);
          int old = maximum_bad_deviation.load(std::memory_order_relaxed);
          while (deviation > old && !maximum_bad_deviation.compare_exchange_weak(
              old, deviation, std::memory_order_relaxed)) {}
          min_atomic(first_bad, index);
        }
        processed.fetch_add(1, std::memory_order_relaxed);
      }
    }
  };
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) pool.emplace_back(worker);
  for (auto& thread : pool) thread.join();
  U bad_index = first_bad.load(std::memory_order_relaxed);
  int bad_frequency = -1, bad_walsh = 0;
  if (bad_index != std::numeric_limits<U>::max()) {
    U pair = bad_index % (n * n);
    int a = static_cast<int>(bad_index / (n * n));
    U i = pair / n, j = pair % n;
    int ignored_deviation = 0;
    check_phase(base ^ extension_mask(seeds[i], seeds[j], a), &bad_frequency, &bad_walsh,
                &ignored_deviation);
  }
  const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();
  std::ofstream output(argv[3]);
  if (!output) return 4;
  output << "{\n"
         << "  \"schema\": \"a28-r016-extension-cpu-v1\",\n"
         << "  \"seed_count\": " << n << ",\n"
         << "  \"threads\": " << threads << ",\n"
         << "  \"phase_count_expected\": " << total << ",\n"
         << "  \"phase_count_processed\": " << processed.load() << ",\n"
         << "  \"full_coverage\": " << (processed.load() == total ? "true" : "false") << ",\n"
         << "  \"f0_six_mask_hex\": \"0x" << std::hex << base << std::dec << "\",\n"
         << "  \"bad_phase_count\": " << bad_phase_count.load() << ",\n"
         << "  \"bad_walsh_term_count\": " << bad_term_count.load() << ",\n"
         << "  \"maximum_bad_abs_deviation_from_8\": " << maximum_bad_deviation.load() << ",\n"
         << "  \"first_bad_index\": " << (bad_index == std::numeric_limits<U>::max() ? -1LL : static_cast<long long>(bad_index)) << ",\n"
         << "  \"first_bad_frequency\": " << bad_frequency << ",\n"
         << "  \"first_bad_walsh\": " << bad_walsh << ",\n"
         << "  \"elapsed_seconds\": " << elapsed << ",\n"
         << "  \"all_abs_walsh_equal_8\": " << (bad_phase_count.load() == 0 && processed.load() == total ? "true" : "false") << "\n"
         << "}\n";
  output.close();
  if (processed.load() != total || bad_phase_count.load() != 0) return 10;
  return 0;
}
