/* Independent full restriction-weight census with an inverted enumeration order. */
#define main flat_census_original_main
#include "flat_census.c"
#undef main

enum { WORDS = X / 64 };
static uint64_t table[SAMPLES][WORDS];
static unsigned pivot_masks[N_PATTERNS];
static uint64_t reverse_prefix[N_PATTERNS + 1];

typedef struct {
    unsigned id, shards;
    uint64_t histogram[SAMPLES][POINTS + 1];
} WeightWorker;

static void generate_truth_by_mobius(const char *input_path) {
    FILE *input = fopen(input_path, "r");
    if (!input) abort();
    for (unsigned sample = 0; sample < SAMPLES; ++sample) {
        char coefficients[TERMS + 2];
        if (fscanf(input, "%131s", coefficients) != 1) abort();
        unsigned values[X] = {0};
        unsigned count = 0, cubic = 0;
        for (unsigned a = 0; a < X; ++a) {
            if (__builtin_popcount(a) > 3) continue;
            assert(coefficients[count] == '0' || coefficients[count] == '1');
            values[a] = coefficients[count++] - '0';
            cubic |= values[a] && __builtin_popcount(a) == 3;
        }
        assert(count == TERMS && cubic && !coefficients[count]);
        for (unsigned bit = 1; bit < X; bit <<= 1) {
            for (unsigned x = 0; x < X; ++x) if (x & bit) values[x] ^= values[x ^ bit];
        }
        for (unsigned x = 0; x < X; ++x) table[sample][x / 64] |= (uint64_t)values[x] << (x % 64);
    }
    fclose(input);
}

static void reversed_rref_index(void) {
    unsigned n = 0;
    for (unsigned pivots = X - 1; pivots > 0; --pivots) {
        if (__builtin_popcount(pivots) != D) continue;
        unsigned freedom = 0;
        for (unsigned c = 0; c < N; ++c) {
            if ((pivots >> c) & 1) continue;
            freedom += __builtin_popcount(pivots & ((1u << c) - 1));
        }
        pivot_masks[n] = pivots;
        reverse_prefix[n + 1] = reverse_prefix[n] + (UINT64_C(1) << freedom);
        ++n;
    }
    assert(n == N_PATTERNS);
    uint64_t gaussian[D + 1] = {1};
    for (unsigned m = 1; m <= N; ++m) {
        for (unsigned k = m < D ? m : D; k > 0; --k) {
            gaussian[k] += (UINT64_C(1) << (m - k)) * gaussian[k - 1];
        }
    }
    assert(gaussian[D] == reverse_prefix[n]);
}

static void *weight_worker(void *opaque) {
    WeightWorker *job = opaque;
    uint64_t total = reverse_prefix[N_PATTERNS];
    uint64_t lo = total * job->id / job->shards, hi = total * (job->id + 1) / job->shards;
    for (unsigned p = 0; p < N_PATTERNS; ++p) {
        unsigned pivot_mask = pivot_masks[p];
        unsigned column_mask = X - 1 - pivot_mask;
        unsigned offsets[FIBRE], n = 0;
        for (unsigned z = 0; z < X; ++z) if (!(z & pivot_mask)) offsets[n++] = z;
        assert(n == FIBRE);
        uint64_t first = lo < reverse_prefix[p] ? reverse_prefix[p] : lo;
        uint64_t last = hi > reverse_prefix[p + 1] ? reverse_prefix[p + 1] : hi;
        for (uint64_t serial = first; serial < last; ++serial) {
            unsigned free_bits = (unsigned)(serial - reverse_prefix[p]);
            unsigned rows[D], r = 0;
            for (unsigned i = 0; i < N; ++i) {
                if (((pivot_mask >> i) & 1) == 0) continue;
                unsigned row = 1u << i;
                for (unsigned j = i + 1; j < N; ++j) {
                    if ((column_mask >> j) & 1) {
                        row |= (free_bits & 1) << j;
                        free_bits >>= 1;
                    }
                }
                rows[r++] = row;
            }
            assert(!free_bits);
            unsigned linear[POINTS] = {0};
            for (unsigned b = 0; b < D; ++b) {
                for (unsigned j = 0; j < (1u << b); ++j) linear[(1u << b) | j] = linear[j] ^ rows[b];
            }
            for (unsigned offset = 0; offset < FIBRE; ++offset) {
                uint64_t flat[WORDS] = {0};
                for (unsigned a = 0; a < POINTS; ++a) {
                    unsigned x = offsets[offset] ^ linear[a];
                    flat[x / 64] |= UINT64_C(1) << (x % 64);
                }
                for (unsigned f = 0; f < SAMPLES; ++f) {
                    unsigned weight = 0;
                    for (unsigned word = 0; word < WORDS; ++word) {
                        weight += __builtin_popcountll(flat[word] & table[f][word]);
                    }
                    ++job->histogram[f][weight];
                }
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc != 9) return 2;
    unsigned shards = (unsigned)strtoul(argv[2], NULL, 10);
    unsigned first = (unsigned)strtoul(argv[3], NULL, 10);
    unsigned threads = (unsigned)strtoul(argv[4], NULL, 10);
    cores = (unsigned)strtoul(argv[5], NULL, 10);
    double planned_start = strtod(argv[6], NULL);
    if (!shards || !threads || first + threads > shards || !cores) return 2;
    generate_truth_by_mobius(argv[1]);
    reversed_rref_index();
    FILE *out = fopen(argv[7], "w");
    usage_file = fopen(argv[8], "w");
    if (!out || !usage_file) abort();
    fputs("epoch,elapsed_seconds,system_active_percent,task_core_equivalents,task_machine_percent\n", usage_file);
    WeightWorker *work = calloc(threads, sizeof(*work));
    pthread_t *ids = calloc(threads, sizeof(*ids));
    if (!work || !ids) abort();
    while (planned_start > now(CLOCK_REALTIME)) sleep_for(planned_start - now(CLOCK_REALTIME));
    double epoch = now(CLOCK_REALTIME), wall = now(CLOCK_MONOTONIC), cpu = now(CLOCK_PROCESS_CPUTIME_ID);
    pthread_t stats_thread;
    if (pthread_create(&stats_thread, NULL, monitor, NULL)) abort();
    for (unsigned i = 0; i < threads; ++i) {
        work[i].id = first + i;
        work[i].shards = shards;
        if (pthread_create(&ids[i], NULL, weight_worker, &work[i])) abort();
    }
    for (unsigned i = 0; i < threads; ++i) if (pthread_join(ids[i], NULL)) abort();
    double used_cpu = now(CLOCK_PROCESS_CPUTIME_ID) - cpu;
    double used_wall = now(CLOCK_MONOTONIC) - wall, end = now(CLOCK_REALTIME);
    atomic_store_explicit(&finished, 1, memory_order_relaxed);
    if (pthread_join(stats_thread, NULL)) abort();

    uint64_t total = reverse_prefix[N_PATTERNS];
    fprintf(out, "{\"start_epoch\":%.6f,\"end_epoch\":%.6f,\"wall_seconds\":%.6f,\"cpu_seconds\":%.6f,"
            "\"logical_cpu_count\":%u,\"worker_count\":%u,\"total_shards\":%u,\"shards\":[",
            epoch, end, used_wall, used_cpu, cores, threads, shards);
    for (unsigned i = 0; i < threads; ++i) {
        unsigned id = work[i].id;
        uint64_t from = total * id / shards, to = total * (id + 1) / shards;
        if (i) fputc(',', out);
        fprintf(out, "{\"id\":%u,\"linear_from\":%" PRIu64 ",\"linear_to\":%" PRIu64
                ",\"affine_count\":%" PRIu64 ",\"histograms\":[", id, from, to, FIBRE * (to - from));
        for (unsigned j = 0; j < SAMPLES; ++j) {
            uint64_t seen = 0;
            fputs(j ? ",[" : "[", out);
            for (unsigned w = 0; w <= POINTS; ++w) {
                uint64_t count = work[i].histogram[j][w];
                if ((w & 1) && count) abort();
                seen += count;
                fprintf(out, "%s%" PRIu64, w ? "," : "", count);
            }
            if (seen != FIBRE * (to - from)) abort();
            fputc(']', out);
        }
        fputs("]}", out);
    }
    fputs("]}\n", out);
    if (fclose(out) || fclose(usage_file)) abort();
    free(ids);
    free(work);
    return 0;
}
