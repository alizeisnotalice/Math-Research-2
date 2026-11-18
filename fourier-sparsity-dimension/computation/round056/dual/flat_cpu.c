#include <inttypes.h>
#include <mach/host_info.h>
#include <mach/mach.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

enum { MAX_WORKERS=64, BINS=193, LOW_LIMIT=4096 };
typedef struct {
    unsigned id;
    unsigned norm;
    unsigned codim;
    unsigned character;
    uint32_t flat;
} Witness;
typedef struct {
    uint64_t hist[BINS];
    uint64_t location_hash;
    unsigned low_count;
    Witness witnesses[LOW_LIMIT];
    double cpu_seconds;
} Totals;
static unsigned first,count,phase_limit;
static uint32_t affine[64];
static FILE *telemetry;
static atomic_uint next_block,active_threads;

static double now(clockid_t id) {
    struct timespec t;
    if(clock_gettime(id,&t))abort();
    return t.tv_sec+t.tv_nsec*1.e-9;
}
static uint32_t truth_of(unsigned id) {
    return ((id&1u)<<3)|((id&14u)<<4)|((id&2032u)<<5)|((id&67106816u)<<6);
}
static uint64_t keyed(unsigned i,unsigned L) {
    uint64_t a=((uint64_t)i<<8)|L;
    a=(a^(a>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    a=(a^(a>>27))*UINT64_C(0x94d049bb133111eb);
    return a^(a>>31);
}
static int flat_rank(uint32_t flat) {
    unsigned anchor=__builtin_ctz(flat);
    unsigned pivot[5]={0}, rank=0;
    for(;flat;flat&=flat-1) {
        unsigned value=(unsigned)__builtin_ctz(flat)^anchor;
        while(value) {
            unsigned index=31u-__builtin_clz(value);
            if(!pivot[index]){pivot[index]=value;rank++;break;}
            value^=pivot[index];
        }
    }
    return (int)rank;
}
static Witness classify(unsigned id,uint32_t mask,unsigned norm) {
    if(mask==0)return (Witness){id,norm,0,0,0};
    for(unsigned b=0;b<64;b++) {
        uint32_t defect=mask^affine[b];
        unsigned size=__builtin_popcount(defect);
        if(!(size>0 && size<=8 && (size&(size-1))==0))continue;
        unsigned k=__builtin_ctz(size);
        if(flat_rank(defect)==(int)k)return (Witness){id,norm,5-k,b,defect};
    }
    return (Witness){id,norm,255,0,0};
}
static unsigned l1(uint32_t mask) {
    int W[32];
    for(unsigned x=0;x<32;x++)W[x]=(mask>>x)&1u ? -1:1;
    for(unsigned h=1;h<32;h*=2) {
        for(unsigned base=0;base<32;base+=2*h)
            for(unsigned x=base;x<base+h;x++){
                int a=W[x],b=W[x+h]; W[x]=a+b; W[x+h]=a-b;
            }
    }
    unsigned L=0;
    for(unsigned x=0;x<32;x++) L+=abs(W[x]);
    return L;
}
static void payload(unsigned id,Totals*t) {
    uint32_t mask=truth_of(id);
    unsigned L=l1(mask);
    if(L>=BINS)abort();
    // In the expanded CPU audit, transform every one of the distinct affine-character variants.
    for(unsigned phase=1;phase<phase_limit;phase++) {
        if(l1(mask^affine[phase])!=L)abort();
    }
    t->hist[L]++;
    t->location_hash^=keyed(id,L);
    if(L<96) {
        if(t->low_count==LOW_LIMIT)abort();
        t->witnesses[t->low_count++]=classify(id,mask,L);
    }
}
static void *worker(void *arg) {
    Totals *t=arg;
    double begin=now(CLOCK_THREAD_CPUTIME_ID);
    for(;;){
        unsigned offset=atomic_fetch_add(&next_block,2048);
        if(offset>=count)break;
        unsigned end=offset+2048;if(end>count)end=count;
        for(unsigned x=offset;x<end;x++)payload(first+x,t);
    }
    t->cpu_seconds=now(CLOCK_THREAD_CPUTIME_ID)-begin;
    atomic_fetch_sub(&active_threads,1);
    return NULL;
}
static host_cpu_load_info_data_t ticks(void) {
    host_cpu_load_info_data_t out;
    mach_msg_type_number_t n=HOST_CPU_LOAD_INFO_COUNT;
    if(host_statistics(mach_host_self(),HOST_CPU_LOAD_INFO,(host_info_t)&out,&n)!=KERN_SUCCESS)abort();
    return out;
}
static void *monitor(void *unused) {
    (void)unused; const struct timespec pause={0,100000000};
    double cpu=now(CLOCK_PROCESS_CPUTIME_ID), wall=now(CLOCK_MONOTONIC);
    host_cpu_load_info_data_t old=ticks();
    while(atomic_load(&active_threads)){
        nanosleep(&pause,NULL);
        double nc=now(CLOCK_PROCESS_CPUTIME_ID),nw=now(CLOCK_MONOTONIC);
        host_cpu_load_info_data_t new=ticks();
        unsigned long all=0,idle=new.cpu_ticks[CPU_STATE_IDLE]-old.cpu_ticks[CPU_STATE_IDLE];
        for(unsigned i=0;i<CPU_STATE_MAX;i++)all+=new.cpu_ticks[i]-old.cpu_ticks[i];
        fprintf(telemetry,"%.6f,%.9f,%.6f,%.6f,%u\n",now(CLOCK_REALTIME),nw-wall,
                all?100.0*(all-idle)/all:-1,100.0*(nc-cpu)/(nw-wall),atomic_load(&active_threads));
        fflush(telemetry);old=new;wall=nw;cpu=nc;
    }
    return NULL;
}
static int compare(const void *a,const void *b) {
    uint32_t x=((const Witness*)a)->id,y=((const Witness*)b)->id;
    return (x>y)-(x<y);
}
int main(int argc,char**argv) {
    if(argc!=8)return 2;
    first=(unsigned)strtoul(argv[1],NULL,10);count=(unsigned)strtoul(argv[2],NULL,10);
    unsigned threads=(unsigned)strtoul(argv[3],NULL,10);
    double planned=strtod(argv[4],NULL);
    phase_limit=(unsigned)strtoul(argv[7],NULL,10);
    if((phase_limit!=1&&phase_limit!=64)||!count||!threads||threads>MAX_WORKERS||(uint64_t)first+count>(1ull<<26))return 2;
    FILE*out=fopen(argv[5],"w");telemetry=fopen(argv[6],"w");
    if(!out||!telemetry)return 3;
    for(unsigned r=0;r<64;r++)
        for(unsigned x=0;x<32;x++)
            affine[r]|=(uint32_t)((__builtin_parity(r&x)^(r>>5))&1)<<x;
    fputs("unix_epoch,sample_seconds,whole_machine_cpu_active_percent,a28_process_cpu_percent,active_workers\n",telemetry);
    while(planned>now(CLOCK_REALTIME)){
        const struct timespec wait={0,5000000};nanosleep(&wait,NULL);
    }
    Totals*per_worker=calloc(threads,sizeof(*per_worker));
    if(!per_worker)abort();
    atomic_init(&active_threads,threads);atomic_init(&next_block,0);
    pthread_t pool[MAX_WORKERS],sampler;
    double start=now(CLOCK_REALTIME), wall=now(CLOCK_MONOTONIC), cpu=now(CLOCK_PROCESS_CPUTIME_ID);
    if(pthread_create(&sampler,NULL,monitor,NULL))abort();
    for(unsigned i=0;i<threads;i++)if(pthread_create(pool+i,NULL,worker,per_worker+i))abort();
    for(unsigned i=0;i<threads;i++)if(pthread_join(pool[i],NULL))abort();
    double end=now(CLOCK_REALTIME),elapsed=now(CLOCK_MONOTONIC)-wall,total_cpu=now(CLOCK_PROCESS_CPUTIME_ID)-cpu;
    if(pthread_join(sampler,NULL))abort();
    uint64_t H[BINS]={0},hash=0,seen=0;
    Witness witnesses[LOW_LIMIT];unsigned used=0;
    for(unsigned i=0;i<threads;i++){
        for(unsigned j=0;j<BINS;j++){H[j]+=per_worker[i].hist[j];seen+=per_worker[i].hist[j];}
        if(used+per_worker[i].low_count>LOW_LIMIT)abort();
        for(unsigned j=0;j<per_worker[i].low_count;j++)witnesses[used++]=per_worker[i].witnesses[j];
        hash^=per_worker[i].location_hash;
    }
    if(seen!=count)abort();
    qsort(witnesses,used,sizeof(*witnesses),compare);
    fprintf(out,"{\"schema\":\"a28-r056-exhaustive-normalized-five-variable-cpu-v1\",\"first_id\":%u,\"count\":%u,"
        "\"threads\":%u,\"start_epoch\":%.6f,\"end_epoch\":%.6f,\"wall\":%.9f,\"cpu_seconds\":%.9f,"
        "\"norm_position_checksum64\":\"%016" PRIx64 "\",\"norm_histogram_units32\":{",
        first,count,threads,start,end,elapsed,total_cpu,hash);
    unsigned printed=0;
    for(unsigned i=0;i<BINS;i++)if(H[i])fprintf(out,"%s\"%u\":%" PRIu64,printed++?",":"",i,H[i]);
    fputs("},\"below_three_records\":[",out);
    for(unsigned i=0;i<used;i++){
        Witness w=witnesses[i];
        fprintf(out,"%s[%u,%u,%u,%u,%u]",i?",":"",w.id,w.norm,w.codim,w.character,w.flat);
    }
    fprintf(out,"],\"affine_variants_transformed_for_each_representative\":%u,"
        "\"distinct_raw_truth_tables_transformed\":%" PRIu64 ",\"worker_seconds\":[",
        phase_limit,(uint64_t)phase_limit*count);
    for(unsigned i=0;i<threads;i++)fprintf(out,"%s%.9f",i?",":"",per_worker[i].cpu_seconds);
    fputs("]}\n",out);
    free(per_worker);
    return fclose(out)||fclose(telemetry);
}
