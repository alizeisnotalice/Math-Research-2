#include "FiberCommon.h"
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
static Word characters[32];

typedef struct {
    uint64_t sum[WIDTH], checksum, dependent_instances, sharp_instances;
    Word max_excess,max_fiber;
    double worker_cpu_seconds;
} Stat;

static uint64_t code_hash(uint64_t x){
    x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);
    return x^(x>>31);
}
static void audit(Word id,Stat *s){
    struct Model m=make_model(id);
    Word V=1u<<m.v,Y=1u<<m.y;
    Word supports[256],costs[256],core[16]={0},row[WIDTH]={0},bad=0;
    int means[256];
    for(Word y=0;y<Y;y++){
        for(Word v=0;v<V;v++){
            struct Local f=local_model(id,m,y,v);
            for(Word u=0;u<16;u++){
                int target=control_mean(m,f,u);
                Word mask=negative_mask(target,f),point=(v<<4)|u,support=0;
                for(Word b=0;b<32;b++)if(bits(mask^characters[b])!=16)support|=1u<<b;
                int dc=32-2*(int)bits(mask);
                if(dc%8)abort();
                means[point]=dc/8;
                costs[point]=private_cost(means[point],support,target,&bad);
                supports[point]=support;
            }
        }
        bad|=inspect_address(m,means,costs,supports,core,row);
    }
    bad|=finish(m,core,row);
    if(bad)abort();
    uint64_t hash=code_hash(id);
    for(Word j=0;j<VALID;j++){
        s->sum[j]+=row[j];
        hash=code_hash(hash^row[j]);
    }
    s->checksum^=hash;
    s->dependent_instances+=(row[EXCESS]>0);
    s->sharp_instances+=row[ACTUAL_SHARP]>0;
    if(row[EXCESS]>s->max_excess)s->max_excess=row[EXCESS];
    if(row[FIBER_MAX]>s->max_fiber)s->max_fiber=row[FIBER_MAX];
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
    for(Word b=0;b<32;b++)for(Word x=0;x<32;x++)characters[b]|=sign_parity(b&x)<<x;
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
        total.dependent_instances+=per[i].dependent_instances;
        total.sharp_instances+=per[i].sharp_instances;
        if(per[i].max_excess>total.max_excess)total.max_excess=per[i].max_excess;
        if(per[i].max_fiber>total.max_fiber)total.max_fiber=per[i].max_fiber;
        for(Word j=0;j<VALID;j++)total.sum[j]+=per[i].sum[j];
    }
    fprintf(out,"{\"first_id\":%u,\"count\":%u,\"threads\":%u,\"start_epoch\":%.6f,\"end_epoch\":%.6f,"
        "\"wall\":%.9f,\"cpu_seconds\":%.9f,\"positional_checksum64\":\"%016" PRIx64 "\",\"dependent_instances\":%" PRIu64
        ",\"sharp_instances\":%" PRIu64 ",\"max_excess\":%u,\"max_fiber\":%u,\"sums\":[",
        first,count,threads,start,end,elapsed,cpu,total.checksum,total.dependent_instances,total.sharp_instances,total.max_excess,total.max_fiber);
    for(Word i=0;i<VALID;i++)fprintf(out,"%s%" PRIu64,i?",":"",total.sum[i]);
    fputs("],\"worker_cpu_seconds\":[",out);
    for(Word i=0;i<threads;i++)fprintf(out,"%s%.9f",i?",":"",per[i].worker_cpu_seconds);
    fputs("]}\n",out);
    return fclose(out)||fclose(telemetry);
}
