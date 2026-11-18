#include <inttypes.h>
#include <mach/host_info.h>
#include <mach/mach.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { MAX_CONFIGS=512, MAX_WORKERS=64, GAP_BINS=256 };
typedef struct { unsigned rank; unsigned rows[5]; unsigned labels[32]; } Config;
typedef struct {
    uint64_t full[7], zero[7], below[7], undersized[7], hist[7][GAP_BINS];
    unsigned example_set[7], examples[7];
    int minimum[7];
    double thread_seconds;
} Totals;
static Config configurations[MAX_CONFIGS];
static unsigned cfg_number;
static uint64_t next_control;
static uint32_t first, count;
static atomic_uint open_workers;
static atomic_uint next_block;
static double epoch_base;
static FILE *telemetry;

static double now(clockid_t input) {
    struct timespec t;
    if(clock_gettime(input,&t)) abort();
    return t.tv_sec+t.tv_nsec*1.e-9;
}
static uint32_t mix(uint32_t x) {
    x ^= x >> 16; x *= UINT32_C(0x7feb352d); x ^= x >> 15; x *= UINT32_C(0x846ca68b);
    return x ^ (x >> 16);
}
static uint32_t translated(uint32_t s,unsigned a) {
    uint32_t out=0;
    while(s) { unsigned i=__builtin_ctz(s); s &= s-1; out |= UINT32_C(1)<<(i^a); }
    return out;
}
static void make_configs(void) {
    uint32_t all[6][200]={{1}};
    unsigned found[6]={1,0,0,0,0,0};
    for(unsigned r=1;r<=5;r++) {
        for(unsigned i=0;i<found[r-1];i++) {
            uint32_t s=all[r-1][i];
            for(unsigned x=1;x<32;x++) {
                if(s&(UINT32_C(1)<<x))continue;
                uint32_t ext=s|translated(s,x);
                unsigned same=0;
                for(unsigned j=0;j<found[r];j++)same|=all[r][j]==ext;
                if(!same)all[r][found[r]++]=ext;
            }
        }
        for(unsigned i=0;i<found[r];i++) {
            Config *c = configurations+cfg_number++;
            c->rank=r; unsigned k=0;
            uint32_t partial=1;
            for(unsigned x=1;x<32;x++) {
                if((all[r][i] & (UINT32_C(1)<<x)) && !(partial&(UINT32_C(1)<<x))) {
                    c->rows[k++]=x; partial |= translated(partial,x);
                }
            }
            if(k!=r)abort();
            for(unsigned x=0;x<32;x++)
                for(unsigned y=0;y<r;y++)c->labels[x]|=(unsigned)__builtin_parity(c->rows[y]&x)<<y;
        }
    }
    if(found[1]!=31||found[2]!=155||found[3]!=155||found[4]!=31||found[5]!=1||cfg_number!=373)abort();
}
static unsigned mass_abs(int x){return x<0 ? (unsigned)(-x):(unsigned)x;}
static void record(Totals*t,unsigned rank,unsigned l1,uint32_t truth) {
    int delta=(int)l1-(64-(64>>rank));
    t->full[rank]++;
    if(delta<0){ t->below[rank]++; return; }
    if(delta>=GAP_BINS)abort();
    t->hist[rank][delta]++;
    if(delta==0){t->zero[rank]++; return;}
    if(delta<(rank>=4?16:8))t->undersized[rank]++;
    if(delta<t->minimum[rank]||(delta==t->minimum[rank]&&truth<t->examples[rank])) {
        t->minimum[rank]=delta; t->examples[rank]=truth; t->example_set[rank]=1;
    }
}
static void audit(uint32_t table,Totals *totals) {
    int W[32];
    for(unsigned x=0;x<32;x++)W[x]=(table>>x)&1 ? -1:1;
    for(unsigned d=1;d<32;d*=2)
        for(unsigned base=0;base<32;base+=2*d)
            for(unsigned x=base;x<base+d;x++){int a=W[x],b=W[x+d];W[x]=a+b;W[x+d]=a-b;}
    unsigned weights[32],sum=0,n=0;
    unsigned nonzero[32];
    for(unsigned x=0;x<32;x++){
        weights[x]=mass_abs(W[x]);sum+=weights[x];
        if(weights[x])nonzero[n++]=x;
    }
    for(unsigned i=0;i<cfg_number;i++) {
        const Config*c=configurations+i;
        unsigned mass[32]={0}; uint32_t active=0;
        for(unsigned x=0;x<n;x++){
            unsigned freq=nonzero[x],label=c->labels[freq];
            mass[label]+=weights[freq]; active|=UINT32_C(1)<<label;
        }
        if(!active)abort();
        // In a private linear projection, a dense set already spans all affine labels.
        const unsigned r=c->rank, m=1u<<r;
        if(__builtin_popcount(active)>m/2) {
            if(r>=2)for(unsigned off=0;off<m;off++)record(totals,r,sum-mass[off],table);
            record(totals,r+1,sum,table);
            continue;
        }
        // Affine hull in the quotient label coordinates.
        unsigned anchor=__builtin_ctz(active);
        uint32_t hull=1;
        for(uint32_t remaining=active;remaining;remaining&=remaining-1) {
            unsigned u=__builtin_ctz(remaining)^anchor;
            if(!(hull&(UINT32_C(1)<<u)))hull|=translated(hull,u);
            if(__builtin_popcount(hull)==m)break;
        }
        unsigned dim=(unsigned)__builtin_popcount(hull);
        if(dim==m) record(totals,r+1,sum,table);
        if(r<2 || dim<m/2)continue;
        if(dim==m) {
            for(unsigned off=0;off<m;off++)record(totals,r,sum-mass[off],table);
        } else {
            uint32_t affine=translated(hull,anchor);
            for(unsigned off=0;off<m;off++)
                if(!(affine&(UINT32_C(1)<<off)))record(totals,r,sum-mass[off],table);
        }
    }
}
static void *work(void *arg) {
    Totals*t=arg;
    double started=now(CLOCK_THREAD_CPUTIME_ID);
    for(unsigned r=0;r<=6;r++)t->minimum[r]=999;
    const unsigned chunk=32;
    for(;;) {
        unsigned from=atomic_fetch_add(&next_block,chunk);
        if(from>=count)break;
        unsigned to=from+chunk; if(to>count)to=count;
        for(unsigned i=from;i<to;i++)audit(mix(first+i),t);
    }
    t->thread_seconds=now(CLOCK_THREAD_CPUTIME_ID)-started;
    atomic_fetch_sub(&open_workers,1);
    return NULL;
}
static host_cpu_load_info_data_t ticks(void) {
    host_cpu_load_info_data_t out;
    mach_msg_type_number_t n=HOST_CPU_LOAD_INFO_COUNT;
    if(host_statistics(mach_host_self(), HOST_CPU_LOAD_INFO,(host_info_t)&out,&n)!=KERN_SUCCESS)abort();
    return out;
}
static void *sample(void *ignored) {
    (void)ignored; const struct timespec pause={0,120000000};
    double cpu=now(CLOCK_PROCESS_CPUTIME_ID),wall=now(CLOCK_MONOTONIC);
    host_cpu_load_info_data_t old=ticks();
    while(atomic_load(&open_workers)){
        nanosleep(&pause,NULL);
        double nc=now(CLOCK_PROCESS_CPUTIME_ID),nw=now(CLOCK_MONOTONIC);
        host_cpu_load_info_data_t new=ticks();
        unsigned long sum=0,idle=new.cpu_ticks[CPU_STATE_IDLE]-old.cpu_ticks[CPU_STATE_IDLE];
        for(unsigned i=0;i<CPU_STATE_MAX;i++)sum+=new.cpu_ticks[i]-old.cpu_ticks[i];
        fprintf(telemetry,"%.6f,%.9f,%.6f,%.6f,%u\n",now(CLOCK_REALTIME),nw-wall,
                sum? 100.0*(sum-idle)/sum : -1,100.0*(nc-cpu)/(nw-wall),atomic_load(&open_workers));
        fflush(telemetry); old=new; wall=nw;cpu=nc;
    }
    return NULL;
}
int main(int argc,char **argv) {
    make_configs();
    if(argc==3 && !strcmp(argv[1],"--export")) {
        FILE*out=fopen(argv[2],"w"); if(!out)return 1;
        fputs("[",out);
        for(unsigned i=0;i<cfg_number;i++){
            Config*c=configurations+i;
            fprintf(out,"%s{\"rank\":%u,\"labels\":[",i?",":"",c->rank);
            for(unsigned j=0;j<32;j++)fprintf(out,"%s%u",j?",":"",c->labels[j]);
            fprintf(out,"]}");
        }
        fputs("]\n",out);return fclose(out)!=0;
    }
    if(argc!=7)return 2;
    first=(uint32_t)strtoul(argv[1],NULL,10); count=(uint32_t)strtoul(argv[2],NULL,10);
    unsigned workers=(unsigned)strtoul(argv[3],NULL,10);
    if(!count||!workers||workers>MAX_WORKERS||(uint64_t)first+count>(1ull<<32))return 2;
    epoch_base=strtod(argv[4],NULL);
    FILE*out=fopen(argv[5],"w"); telemetry=fopen(argv[6],"w");
    if(!out||!telemetry) return 3;
    fprintf(telemetry,"unix_epoch,sample_seconds,whole_machine_cpu_active_percent,a28_process_cpu_percent,active_workers\n");
    while(epoch_base>now(CLOCK_REALTIME)) {
        const struct timespec wait={0,5000000}; nanosleep(&wait,NULL);
    }
    double start=now(CLOCK_REALTIME), wall=now(CLOCK_MONOTONIC), cpu=now(CLOCK_PROCESS_CPUTIME_ID);
    pthread_t worker[MAX_WORKERS], sampler; Totals all[MAX_WORKERS]={0};
    atomic_init(&open_workers,workers); atomic_init(&next_block,0);
    if(pthread_create(&sampler,NULL,sample,NULL))abort();
    for(unsigned i=0;i<workers;i++)if(pthread_create(worker+i,NULL,work,all+i))abort();
    for(unsigned i=0;i<workers;i++)if(pthread_join(worker[i],NULL))abort();
    double end=now(CLOCK_REALTIME),elapsed=now(CLOCK_MONOTONIC)-wall,cpu_used=now(CLOCK_PROCESS_CPUTIME_ID)-cpu;
    if(pthread_join(sampler,NULL))abort();
    fprintf(out,"{\"schema\":\"a28-r055-five-variable-full-affine-evaluation-cpu-v1\","
      "\"first_id\":%u,\"count\":%u,\"threads\":%u,\"configurations\":%u,\"start_epoch\":%.6f,"
      "\"end_epoch\":%.6f,\"wall\":%.8f,\"cpu_seconds\":%.8f,\"by_rank\":{",
      first,count,workers,cfg_number,start,end,elapsed,cpu_used);
    for(unsigned r=2;r<=6;r++){
        uint64_t full=0,zero=0,below=0,bad=0,hist[GAP_BINS]={0};
        uint32_t example=0; int minimum=999;
        for(unsigned j=0;j<workers;j++){
            Totals*t=all+j;full+=t->full[r];zero+=t->zero[r];below+=t->below[r];bad+=t->undersized[r];
            for(unsigned k=0;k<GAP_BINS;k++)hist[k]+=t->hist[r][k];
            if(t->example_set[r]&&(t->minimum[r]<minimum||(t->minimum[r]==minimum&&t->examples[r]<example))){
                example=t->examples[r];minimum=t->minimum[r];
            }
        }
        fprintf(out,"%s\"%u\":{\"full\":%" PRIu64 ",\"endpoint\":%" PRIu64 ",\"below_bound\":%" PRIu64
          ",\"positive_below_requested_gap\":%" PRIu64 ",\"minimum_positive_units32\":%d,\"minimum_truth_mask\":%u,\"gap_histogram_units32\":{",
            r>2?",":"",r,full,zero,below,bad,minimum==999?0:minimum,example);
        unsigned wrote=0;
        for(unsigned k=0;k<GAP_BINS;k++)if(hist[k])fprintf(out,"%s\"%u\":%" PRIu64,wrote++?",":"",k,hist[k]);
        fputs("}}",out);
    }
    fputs("},\"worker_cpu_seconds\":[",out);
    for(unsigned i=0;i<workers;i++)fprintf(out,"%s%.8f",i?",":"",all[i].thread_seconds);
    fputs("]}\n",out);
    return fclose(out)||fclose(telemetry);
}
