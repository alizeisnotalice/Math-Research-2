#define _DARWIN_C_SOURCE
#include <assert.h>
#include <inttypes.h>
#include <mach/mach.h>
#include <mach/mach_host.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

enum { N = 9, D = 5, X = 1 << N, FIBRE = 1 << (N - D), POINTS = 1 << D,
       N_PATTERNS = 126, SAMPLES = 64, TERMS = 130 };

typedef struct {
    uint64_t start, limit;
    unsigned rows[D][FIBRE];
    unsigned shift[D], mask[D];
    unsigned cosets[FIBRE];
} Pattern;

typedef struct {
    uint64_t witness_id;
    unsigned rows[D];
    unsigned offset, sign_bit;
} Witness;

typedef struct {
    unsigned id, shard_count;
    uint64_t from, to, flat_count;
    uint64_t counts[SAMPLES];
    Witness witnesses[SAMPLES];
} Worker;

static uint64_t truth[X];
static Pattern patterns[N_PATTERNS];
static uint64_t n_subspaces;
static unsigned ctz[POINTS];
static atomic_bool finished;
static FILE *usage_file;
static unsigned cores;

static double now(clockid_t clock) {
    struct timespec t;
    if (clock_gettime(clock, &t)) abort();
    return (double)t.tv_sec + t.tv_nsec * 1e-9;
}

static void sleep_for(double seconds) {
    if (seconds <= 0) return;
    time_t integral = (time_t)seconds;
    struct timespec t = {integral, (long)((seconds - integral) * 1e9)};
    nanosleep(&t, NULL);
}

static int query_ticks(uint32_t ticks[CPU_STATE_MAX]) {
    host_cpu_load_info_data_t cpu;
    mach_msg_type_number_t n = HOST_CPU_LOAD_INFO_COUNT;
    host_t host = mach_host_self();
    kern_return_t rc = host_statistics(host, HOST_CPU_LOAD_INFO, (host_info_t)&cpu, &n);
    mach_port_deallocate(mach_task_self(), host);
    if (rc != KERN_SUCCESS) return 0;
    for (unsigned s = 0; s < CPU_STATE_MAX; ++s) ticks[s] = cpu.cpu_ticks[s];
    return 1;
}

static void *monitor(void *argument) {
    (void)argument;
    uint32_t old[CPU_STATE_MAX], current[CPU_STATE_MAX];
    if (!query_ticks(old)) return NULL;
    double wall = now(CLOCK_MONOTONIC), cpu = now(CLOCK_PROCESS_CPUTIME_ID), start = wall;
    while (!atomic_load_explicit(&finished, memory_order_relaxed)) {
        sleep_for(0.04);
        if (!query_ticks(current)) return NULL;
        double next_wall = now(CLOCK_MONOTONIC), next_cpu = now(CLOCK_PROCESS_CPUTIME_ID);
        uint64_t total = 0;
        for (unsigned s = 0; s < CPU_STATE_MAX; ++s) total += (uint32_t)(current[s] - old[s]);
        uint32_t idle = current[CPU_STATE_IDLE] - old[CPU_STATE_IDLE];
        if (total && next_wall > wall) {
            double task_cores = (next_cpu - cpu) / (next_wall - wall);
            fprintf(usage_file, "%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    now(CLOCK_REALTIME), next_wall - start, 100.0 * (1.0 - (double)idle / total),
                    task_cores, 100.0 * task_cores / cores);
            fflush(usage_file);
        }
        for (unsigned s = 0; s < CPU_STATE_MAX; ++s) old[s] = current[s];
        wall = next_wall;
        cpu = next_cpu;
    }
    return NULL;
}

static void load_samples(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) abort();
    for (unsigned i = 0; i < SAMPLES; ++i) {
        char bits[TERMS + 2];
        if (fscanf(f, "%131s", bits) != 1) abort();
        unsigned monomials[TERMS], count = 0, term = 0, degree3 = 0;
        for (unsigned m = 0; m < X; ++m) {
            unsigned degree = __builtin_popcount(m);
            if (degree > 3) continue;
            if (bits[term] != '0' && bits[term] != '1') abort();
            if (bits[term++] == '1') {
                monomials[count++] = m;
                degree3 |= degree == 3;
            }
        }
        if (term != TERMS || bits[term] != '\0' || !degree3) abort();
        for (unsigned x = 0; x < X; ++x) {
            unsigned bit = 0;
            for (unsigned j = 0; j < count; ++j) bit ^= (x & monomials[j]) == monomials[j];
            truth[x] |= (uint64_t)bit << i;
        }
    }
    char extra[2];
    if (fscanf(f, "%1s", extra) != EOF) abort();
    fclose(f);
}

static void init_patterns(void) {
    unsigned next = 0;
    for (unsigned pivot_mask = 0; pivot_mask < X; ++pivot_mask) {
        if (__builtin_popcount(pivot_mask) != D) continue;
        Pattern *p = &patterns[next++];
        unsigned nonpivots = X - 1 - pivot_mask;
        unsigned entry = 0, row = 0;
        for (unsigned j = 0; j < N; ++j) {
            if (!(pivot_mask & (1u << j))) continue;
            unsigned allowed = nonpivots & ~((1u << (j + 1)) - 1u);
            unsigned n = 0;
            p->shift[row] = entry;
            for (unsigned subset = allowed;; subset = (subset - 1) & allowed) {
                p->rows[row][n++] = (1u << j) | subset;
                if (!subset) break;
            }
            p->mask[row++] = n - 1;
            entry += __builtin_ctz(n);
        }
        p->start = n_subspaces;
        p->limit = n_subspaces += (uint64_t)1 << entry;
        unsigned off = 0;
        for (unsigned subset = nonpivots;; subset = (subset - 1) & nonpivots) {
            p->cosets[off++] = subset;
            if (!subset) break;
        }
        assert(row == D && off == FIBRE);
    }
    assert(next == N_PATTERNS);
    for (unsigned i = 1; i < POINTS; ++i) ctz[i] = __builtin_ctz(i);
    uint64_t numerator = 1, denominator = 1;
    for (unsigned i = 0; i < D; ++i) {
        numerator *= (1u << N) - (1u << i);
        denominator *= (1u << D) - (1u << i);
    }
    assert(numerator % denominator == 0 && numerator / denominator == n_subspaces);
}

static void *census(void *arg) {
    Worker *w = arg;
    w->from = n_subspaces * w->id / w->shard_count;
    w->to = n_subspaces * (w->id + 1) / w->shard_count;
    for (unsigned i = 0; i < SAMPLES; ++i) w->witnesses[i].witness_id = UINT64_MAX;
    for (unsigned i = 0; i < N_PATTERNS; ++i) {
        const Pattern *p = &patterns[i];
        uint64_t lo = w->from > p->start ? w->from : p->start;
        uint64_t hi = w->to < p->limit ? w->to : p->limit;
        for (uint64_t r = lo; r < hi; ++r) {
            unsigned rows[D], flat[POINTS] = {0};
            unsigned within = (unsigned)(r - p->start);
            for (unsigned j = 0; j < D; ++j) {
                rows[j] = p->rows[j][(within >> p->shift[j]) & p->mask[j]];
            }
            for (unsigned a = 1; a < POINTS; ++a) flat[a] = flat[a & (a - 1)] ^ rows[ctz[a]];
            for (unsigned j = 0; j < FIBRE; ++j) {
                unsigned offset = p->cosets[j];
                uint64_t all_ones = UINT64_MAX, some_ones = 0;
                for (unsigned a = 0; a < POINTS; ++a) {
                    uint64_t value = truth[flat[a] ^ offset];
                    all_ones &= value;
                    some_ones |= value;
                }
                uint64_t monochromatic = all_ones | ~some_ones;
                while (monochromatic) {
                    unsigned sample = __builtin_ctzll(monochromatic);
                    monochromatic &= monochromatic - 1;
                    Witness *hit = &w->witnesses[sample];
                    ++w->counts[sample];
                    if (hit->witness_id == UINT64_MAX) {
                        hit->witness_id = FIBRE * r + j;
                        hit->offset = offset;
                        hit->sign_bit = (truth[offset] >> sample) & 1;
                        for (unsigned b = 0; b < D; ++b) hit->rows[b] = rows[b];
                    }
                }
            }
        }
    }
    w->flat_count = FIBRE * (w->to - w->from);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc != 9) return 2;
    unsigned shards = (unsigned)strtoul(argv[2], NULL, 10);
    unsigned first = (unsigned)strtoul(argv[3], NULL, 10);
    unsigned threads = (unsigned)strtoul(argv[4], NULL, 10);
    cores = (unsigned)strtoul(argv[5], NULL, 10);
    double start_time = strtod(argv[6], NULL);
    if (!shards || !threads || first + threads > shards || !cores) return 2;
    load_samples(argv[1]);
    init_patterns();
    FILE *out = fopen(argv[7], "w");
    usage_file = fopen(argv[8], "w");
    if (!out || !usage_file) abort();
    fprintf(usage_file, "epoch,elapsed_seconds,system_active_percent,task_core_equivalents,task_machine_percent\n");
    Worker *workers = calloc(threads, sizeof(*workers));
    pthread_t *threads_ids = calloc(threads, sizeof(*threads_ids));
    if (!workers || !threads_ids) abort();
    while (start_time > now(CLOCK_REALTIME)) sleep_for(start_time - now(CLOCK_REALTIME));
    double epoch = now(CLOCK_REALTIME), wall = now(CLOCK_MONOTONIC), cpu = now(CLOCK_PROCESS_CPUTIME_ID);
    pthread_t stats_thread;
    if (pthread_create(&stats_thread, NULL, monitor, NULL)) abort();
    for (unsigned i = 0; i < threads; ++i) {
        workers[i].id = first + i;
        workers[i].shard_count = shards;
        if (pthread_create(&threads_ids[i], NULL, census, &workers[i])) abort();
    }
    for (unsigned i = 0; i < threads; ++i) if (pthread_join(threads_ids[i], NULL)) abort();
    double used_cpu = now(CLOCK_PROCESS_CPUTIME_ID) - cpu;
    double used_wall = now(CLOCK_MONOTONIC) - wall;
    double end = now(CLOCK_REALTIME);
    atomic_store_explicit(&finished, 1, memory_order_relaxed);
    if (pthread_join(stats_thread, NULL)) abort();

    fprintf(out, "{\"worker_count\":%u,\"logical_cpu_count\":%u,\"total_shards\":%u,"
            "\"first_shard\":%u,\"start_epoch\":%.6f,\"end_epoch\":%.6f,"
            "\"wall_seconds\":%.6f,\"cpu_seconds\":%.6f,\"linear_count\":%" PRIu64
            ",\"full_count\":%" PRIu64 ",\"shards\":[",
            threads, cores, shards, first, epoch, end, used_wall, used_cpu,
            n_subspaces, n_subspaces * FIBRE);
    for (unsigned i = 0; i < threads; ++i) {
        const Worker *w = &workers[i];
        if (i) fputc(',', out);
        fprintf(out, "{\"id\":%u,\"linear_from\":%" PRIu64 ",\"linear_to\":%" PRIu64
                ",\"affine_count\":%" PRIu64 ",\"samples\":[", w->id, w->from, w->to, w->flat_count);
        for (unsigned j = 0; j < SAMPLES; ++j) {
            const Witness *hit = &w->witnesses[j];
            if (j) fputc(',', out);
            fprintf(out, "{\"count\":%" PRIu64 ",\"witness\":", w->counts[j]);
            if (!w->counts[j]) {
                fputs("null}", out);
                continue;
            }
            fprintf(out, "{\"flat_id\":%" PRIu64 ",\"offset\":%u,\"bit\":%u,\"rows\":[",
                    hit->witness_id, hit->offset, hit->sign_bit);
            for (unsigned b = 0; b < D; ++b) fprintf(out, "%s%u", b ? "," : "", hit->rows[b]);
            fputs("]}}", out);
        }
        fputs("]}", out);
    }
    fputs("]}\n", out);
    if (fclose(out) || fclose(usage_file)) abort();
    free(workers);
    free(threads_ids);
    return 0;
}
