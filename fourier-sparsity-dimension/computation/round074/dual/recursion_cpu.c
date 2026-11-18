#include "RecursionCommon.h"
#include <inttypes.h>
#include <mach/mach.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

enum { MAX_WORKERS=64 };
static Word first,count;
static atomic_uint next_block,workers;
static FILE *telemetry;
static Wide characters[64];

typedef struct {
    uint64_t sum[WIDTH], checksum;
    double worker_cpu_seconds;
} Stat;

static uint64_t code_hash(uint64_t x){
    x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);
    return x^(x>>31);
}
static void audit(Word id,Stat *s){
    struct Model m=make_model(id);
    Word B=1u<<m.k,Z=1u<<m.n,row[WIDTH]={0};
    int child[256],parent[256];
    Wide restriction=m.n==6?UINT64_MAX:UINT32_MAX;
    for(Word b=0;b<B;b++){
        struct Leaf f=leaf(id,m,b);
        Wide negative=0;
        for(Word x=0;x<Z;x++){
            Word v=evaluate(f,x);
            negative|=(Wide)v<<x;
            parent[b*Z+x]=1-2*(int)v;
        }
        for(Word beta=0;beta<Z;beta++){
            int sum=(int)__builtin_popcountll((negative^characters[beta])&restriction);
            child[b*Z+beta]=(int)Z-2*sum;
        }
    }
    transform(parent,B*Z);
    inspect(m,child,parent,row);
    if(!row[VALID])abort();
    uint64_t hash=code_hash(id);
    for(Word j=0;j<VALID;j++){
        s->sum[j]+=row[j];
        hash=code_hash(hash^row[j]);
    }
    s->checksum^=hash;
}
static double now(clockid_t id){struct timespec t;if(clock_gettime(id,&t))abort();return t.tv_sec+t.tv_nsec*1e-9;}
static void *work(void *arg){
    Stat*s=arg;double start=now(CLOCK_THREAD_CPUTIME_ID);
    for(;;){
        Word from=atomic_fetch_add(&next_block,32);
        if(from>=count)break;
        Word end=from+32;if(end>count)end=count;
        for(Word i=from;i<end;i++)audit(first+i,s);
    }
    s->worker_cpu_seconds=now(CLOCK_THREAD_CPUTIME_ID)-start;
    atomic_fetch_sub(&workers,1);
    return NULL;
}
static host_cpu_load_info_data_t ticks(void){
    host_cpu_load_info_data_t sample;
    mach_msg_type_number_t length=HOST_CPU_LOAD_INFO_COUNT;
    if(host_statistics(mach_host_self(),HOST_CPU_LOAD_INFO,(host_info_t)&sample,&length)!=KERN_SUCCESS)abort();
    return sample;
}
static void *monitor(void *unused){
    (void)unused;
    const struct timespec delay={0,120000000};
    double cpu=now(CLOCK_PROCESS_CPUTIME_ID),wall=now(CLOCK_MONOTONIC);
    host_cpu_load_info_data_t old=ticks();
    while(atomic_load(&workers)){
        nanosleep(&delay,NULL);
        double new_cpu=now(CLOCK_PROCESS_CPUTIME_ID),new_wall=now(CLOCK_MONOTONIC);
        host_cpu_load_info_data_t fresh=ticks();
        unsigned long all=0,idle=fresh.cpu_ticks[CPU_STATE_IDLE]-old.cpu_ticks[CPU_STATE_IDLE];
        for(Word i=0;i<CPU_STATE_MAX;i++)all+=fresh.cpu_ticks[i]-old.cpu_ticks[i];
        fprintf(telemetry,"%.6f,%.9f,%.6f,%.6f,%u\n",now(CLOCK_REALTIME),new_wall-wall,
            all?100.0*(all-idle)/all:-1,100.0*(new_cpu-cpu)/(new_wall-wall),atomic_load(&workers));
        fflush(telemetry);old=fresh;wall=new_wall;cpu=new_cpu;
    }
    return NULL;
}
int main(int argc,char **argv){
    if(argc!=7)return 2;
    first=(Word)strtoul(argv[1],NULL,10);count=(Word)strtoul(argv[2],NULL,10);
    Word threads=(Word)strtoul(argv[3],NULL,10);double planned=strtod(argv[4],NULL);
    if(!count||!threads||threads>MAX_WORKERS)return 2;
    for(Word b=0;b<64;b++)for(Word x=0;x<64;x++)characters[b]|=(Wide)parity(b&x)<<x;
    FILE*out=fopen(argv[5],"w");telemetry=fopen(argv[6],"w");if(!out||!telemetry)return 3;
    fputs("unix_epoch,sample_seconds,whole_machine_cpu_active_percent,a28_process_cpu_percent,active_workers\n",telemetry);
    while(planned>now(CLOCK_REALTIME)){const struct timespec delay={0,5000000};nanosleep(&delay,NULL);}
    atomic_init(&workers,threads);atomic_init(&next_block,0);
    pthread_t work_ids[MAX_WORKERS],sampler;Stat per[MAX_WORKERS]={0},total={0};
    double elapsed=now(CLOCK_MONOTONIC),start=now(CLOCK_REALTIME),cpu=now(CLOCK_PROCESS_CPUTIME_ID);
    if(pthread_create(&sampler,NULL,monitor,NULL))abort();
    for(Word i=0;i<threads;i++)if(pthread_create(&work_ids[i],NULL,work,per+i))abort();
    for(Word i=0;i<threads;i++)if(pthread_join(work_ids[i],NULL))abort();
    double end=now(CLOCK_REALTIME);elapsed=now(CLOCK_MONOTONIC)-elapsed;cpu=now(CLOCK_PROCESS_CPUTIME_ID)-cpu;
    if(pthread_join(sampler,NULL))abort();
    for(Word i=0;i<threads;i++){
        total.checksum^=per[i].checksum;
        for(Word j=0;j<VALID;j++)total.sum[j]+=per[i].sum[j];
    }
    fprintf(out,"{\"first_id\":%u,\"count\":%u,\"threads\":%u,\"start_epoch\":%.6f,\"end_epoch\":%.6f,"
        "\"wall\":%.9f,\"cpu_seconds\":%.9f,\"positional_checksum64\":\"%016" PRIx64 "\",\"sums\":[",
        first,count,threads,start,end,elapsed,cpu,total.checksum);
    for(Word i=0;i<VALID;i++)fprintf(out,"%s%" PRIu64,i?",":"",total.sum[i]);
    fputs("],\"worker_cpu_seconds\":[",out);
    for(Word i=0;i<threads;i++)fprintf(out,"%s%.9f",i?",":"",per[i].worker_cpu_seconds);
    fputs("]}\n",out);
    return fclose(out)||fclose(telemetry);
}
