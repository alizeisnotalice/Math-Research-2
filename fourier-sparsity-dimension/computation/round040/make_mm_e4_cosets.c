#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define CAPACITY 184320
#define HASH_CAPACITY (1u << 19)
#define MASK UINT64_C(0x1111111111111111)
#define E4_0 UINT64_C(0x0101010101010101)
#define E4_1 UINT64_C(0x1010101010101010)
#define E4_2 UINT64_C(0x1111111100000000)

enum { EXTERIOR_LIMIT = 24, THREAD_LIMIT = 64 };

struct candidate {
    uint64_t rep;
    uint64_t outside;
};

struct thread_result {
    unsigned thread;
    uint64_t hist[25];
    uint64_t first[25], second[25];
    uint64_t pairs_checked;
    uint64_t pairs_passed_weight;
    double cpu_seconds;
};

static struct candidate candidates[CAPACITY];
static uint64_t table[HASH_CAPACITY];
static size_t num_candidates = 0;
static unsigned valid_permutations = 0;
static atomic_size_t next_row;
static uint64_t alternates[8];
static unsigned histogram[49];

static uint64_t quotient(uint64_t h) {
    if (h & (UINT64_C(1)<<60)) h ^= E4_1;
    if (h & (UINT64_C(1)<<56)) h ^= E4_0;
    if (h & (UINT64_C(1)<<28)) h ^= E4_2 ^ E4_0 ^ E4_1;
    return h;
}

static unsigned hash_start(uint64_t x) {
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return (unsigned)x & (HASH_CAPACITY - 1);
}

static int contains_or_add(uint64_t x, int add) {
    unsigned i = hash_start(x);
    while (table[i] != 0 && table[i] != x) i = (i+1) & (HASH_CAPACITY-1);
    if (table[i] == x) return 1;
    if (add) table[i] = x;
    return 0;
}

static int compare(const void *a, const void *b) {
    const struct candidate *x = a, *y = b;
    return (x->rep > y->rep) - (x->rep < y->rep);
}

static void handle_permutation(const unsigned pi[8]) {
    if ((pi[0] ^ pi[4]) != 4) return;
    ++valid_permutations;
    uint64_t slope = 0;
    for (unsigned y=0; y<8; ++y) {
        for (unsigned x=0; x<8; ++x) {
            if (__builtin_parity((pi[x] ^ x) & y))
                slope |= UINT64_C(1) << (x + 8*y);
        }
    }
    for (unsigned c=0; c<256; ++c) {
        uint64_t constant = c * UINT64_C(0x0101010101010101);
        uint64_t h = slope ^ constant;
        uint64_t outside = h & ~MASK;
        const unsigned wt = __builtin_popcountll(outside);
        ++histogram[wt];
        if (wt > EXTERIOR_LIMIT) continue;
        uint64_t rep = quotient(h);
        uint64_t expected = h;
        for (unsigned i=0; i<8; ++i) if ((h^alternates[i]) < expected) expected = h^alternates[i];
        if (rep != expected) abort();
        if (rep == 0 || contains_or_add(rep, 1)) continue;
        if (num_candidates == CAPACITY) abort();
        candidates[num_candidates++] = (struct candidate){rep, outside};
    }
}

static void permutations(unsigned pi[8], unsigned index, unsigned used) {
    if (index == 8) {
        handle_permutation(pi);
        return;
    }
    for (unsigned value=0; value<8; ++value) {
        if (used & (1u<<value)) continue;
        pi[index] = value;
        permutations(pi, index+1, used | (1u<<value));
    }
}

static double clock_seconds(clockid_t clock_id) {
    struct timespec value;
    if (clock_gettime(clock_id, &value) != 0) abort();
    return value.tv_sec + value.tv_nsec / 1e9;
}

static void *search(void *argument) {
    struct thread_result *result = argument;
    const double started = clock_seconds(CLOCK_THREAD_CPUTIME_ID);
    const size_t block = 16;
    for (;;) {
        const size_t start = atomic_fetch_add(&next_row, block);
        if (start >= num_candidates) break;
        size_t stop = start + block;
        if (stop > num_candidates) stop = num_candidates;
        for (size_t i=start; i<stop; ++i) {
            struct candidate left = candidates[i];
            result->pairs_checked += num_candidates - i - 1;
            for (size_t j=i+1; j<num_candidates; ++j) {
                struct candidate right = candidates[j];
                uint64_t other = left.rep ^ right.rep;
                if (other <= right.rep) continue;
                unsigned union_weight = __builtin_popcountll(left.outside | right.outside);
                if (union_weight > EXTERIOR_LIMIT) continue;
                ++result->pairs_passed_weight;
                if (!contains_or_add(other, 0)) continue;
                ++result->hist[union_weight];
                if (!result->first[union_weight]) {
                    result->first[union_weight] = left.rep;
                    result->second[union_weight] = right.rep;
                }
            }
        }
    }
    result->cpu_seconds = clock_seconds(CLOCK_THREAD_CPUTIME_ID) - started;
    return NULL;
}

int main(int argc, char **argv) {
    if (argc != 2) return 1;
    const unsigned count = (unsigned)strtoul(argv[1], NULL, 10);
    if (!count || count > THREAD_LIMIT) return 1;
    for (unsigned i=0; i<8; ++i) {
        if (i&1) alternates[i] ^= E4_0;
        if (i&2) alternates[i] ^= E4_1;
        if (i&4) alternates[i] ^= E4_2;
    }
    unsigned pi[8];
    permutations(pi, 0, 0);
    qsort(candidates, num_candidates, sizeof(struct candidate), compare);
    if (valid_permutations != 5760 || histogram[8] != 720 || histogram[0] != 8) abort();
    fprintf(stderr, "Starting %u worker threads over %zu nonzero cosets\n", count, num_candidates);

    atomic_init(&next_row, 0);
    pthread_t threads[THREAD_LIMIT];
    struct thread_result results[THREAD_LIMIT] = {0};
    const double started = clock_seconds(CLOCK_MONOTONIC);
    for (unsigned i=0; i<count; ++i) {
        results[i].thread = i;
        if (pthread_create(&threads[i], NULL, search, &results[i]) != 0) abort();
    }
    for (unsigned i=0; i<count; ++i) if (pthread_join(threads[i], NULL) != 0) abort();
    double elapsed = clock_seconds(CLOCK_MONOTONIC) - started;
    printf("{\n\"schema\":\"a28-r037-mm-five-dim-e4-v1\",\n");
    printf("\"expected_reproduction\":\"./r37_mm_five <logical_core_count>\",\n");
    printf("\"worker_threads\":%u,\"exterior_weight_limit\":%u,\"compatible_permutations\":%u,\n",
           count, EXTERIOR_LIMIT, valid_permutations);
    printf("\"nonzero_coset_count\":%zu, \"search_elapsed_seconds\":%.9f,\n", num_candidates, elapsed);
    uint64_t pairs = 0, passed = 0;
    double cpu_seconds = 0;
    printf("\"worker_cpu_seconds\":[");
    for (unsigned i=0; i<count; ++i) {
        pairs += results[i].pairs_checked;
        passed += results[i].pairs_passed_weight;
        cpu_seconds += results[i].cpu_seconds;
        printf("%s%.9f", i ? "," : "", results[i].cpu_seconds);
    }
    if (pairs != ((uint64_t)num_candidates * (num_candidates-1))/2) abort();
    printf("],\n\"search_cpu_seconds\":%.9f,\"pairs_checked\":%" PRIu64 ","
           "\"canonical_pairs_passing_exterior_weight\":%" PRIu64 ",\n\"extensions\":[",
           cpu_seconds, pairs, passed);
    int wrote = 0;
    for (unsigned wt=0; wt<=EXTERIOR_LIMIT; ++wt) {
        uint64_t total = 0, a = UINT64_MAX, b = UINT64_MAX;
        for (unsigned i=0; i<count; ++i) {
            total += results[i].hist[wt];
            if (results[i].first[wt] &&
                (results[i].first[wt] < a ||
                 (results[i].first[wt] == a && results[i].second[wt] < b))) {
                a = results[i].first[wt];
                b = results[i].second[wt];
            }
        }
        if (total == 0) continue;
        printf("%s{\"strict_exterior_union_size\":%u,\"joint_support\":%u,"
               "\"number_of_5d_extensions\":%" PRIu64 ",\"witness\":[\"%016" PRIx64 "\",\"%016" PRIx64 "\"]}",
               wrote++ ? "," : "", wt, wt+16, total, a, b);
    }
    printf("]\n}\n");
    return 0;
}
