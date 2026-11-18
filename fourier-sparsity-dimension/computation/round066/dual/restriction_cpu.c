#include "RestrictionCommon.h"
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

typedef struct {
    uint64_t sum[WIDTH], checksum, restricted_rank_losses, collisions, missed, collision_by_hk[2][5];
    uint64_t tested_by_hk[2][5];
    double worker_cpu_seconds;
} Stat;

static void transform(int *a,Word length){
    for(Word step=1;step<length;step<<=1)
        for(Word p=0;p<length;p+=2*step)
            for(Word j=0;j<step;j++){
                int u=a[p+j],v=a[p+j+step];
                a[p+j]=u+v; a[p+j+step]=u-v;
            }
}
static uint64_t code_hash(uint64_t x){
    x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);
    return x^(x>>31);
}
static void audit(Word id,Stat *s){
    struct Model m=create_model(id);
    Word N=1u<<m.n,M=1u<<m.h,V=1u<<m.k,small=8*V,scale=M/V;
    int full[512],actual[256],folded[256]={0};
    Word result[WIDTH]={0};
    result[H]=m.h;
    result[K]=m.k;
    result[HOLED]=m.holes;
    for(Word x=0;x<N;x++)full[x]=evaluate(x,m);
    for(Word x=0;x<small;x++)actual[x]=evaluate(restriction_point(x,m),m);
    transform(full,N);
    transform(actual,small);
    Word basis_before[3]={0},basis_after[3]={0};
    Word step=1u<<(m.n-m.h+1+m.holes),spectral_hash=0;
    for(Word beta=0;beta<8;beta++){
        Word first_tail=M,band_size=0,visible=0;
        for(Word a=0;a<M;a++){
            int c=full[beta*M+a];
            int abs_c=abs(c);
            Word pattern=core_pattern_weight(a,m);
            Word expected_odd=beta==0 && (!m.holes || (pattern&1u));
            if(abs_c%step || ((abs_c/step)&1u)!=expected_odd)abort();
            result[ODD_CORE]+=expected_odd;
            if(!beta){
                if((c!=0)!=(pattern!=2))abort();
                result[CORE]+=(c!=0);
            }
            folded[beta*V+image(a,m)]+=parity(a&m.shift)?-c:c;
            if(!c)continue;
            result[ORIG]++;
            band_size++;
            if(first_tail==M)first_tail=a;
        }
        for(Word a=0;a<V;a++){
            int c=actual[beta*V+a];
            if(folded[beta*V+a]!=c*(int)scale)abort();
            spectral_hash=spectral_hash*16777619u^(Word)c;
            if(c){result[RESTRICTED]++;visible=1;}
        }
        if(!beta)continue;
        if(visible)result[D_AFTER]+=add_rank(beta,basis_after,3);
        if(!band_size||!add_rank(beta,basis_before,3))continue;
        result[D_BEFORE]++;
        result[MISSING]+=!visible;
        Word proj=image(first_tail,m),pairs=0;
        result[ELIGIBLE]+=band_size-1;
        for(Word a=first_tail+1;a<M;a++)if(full[beta*M+a]&&image(a,m)==proj)pairs++;
        result[PAIRS]+=pairs;
        result[COLLIDING]+=(pairs>0);
        if(!visible&&!pairs)abort();
    }
    if(result[CORE]!=(m.holes?5*M/8:M) || result[ODD_CORE]!=(m.holes?M/2:M) || result[D_AFTER]>result[D_BEFORE] ||
       result[D_BEFORE]-result[D_AFTER]>result[MISSING] ||
       result[MISSING]>result[COLLIDING] || result[COLLIDING]>result[PAIRS] ||
       result[PAIRS]>result[ELIGIBLE] || result[ELIGIBLE]>result[ORIG]-result[CORE]-result[D_BEFORE])abort();
    result[SPECTRAL_HASH]=spectral_hash;
    uint64_t hash=code_hash(id);
    for(Word j=0;j<VALID;j++){
        s->sum[j]+=result[j];
        hash=code_hash(hash^result[j]);
    }
    s->checksum^=hash;
    s->restricted_rank_losses+=(result[D_BEFORE]!=result[D_AFTER]);
    s->collisions+=(result[COLLIDING]>0);
    s->missed+=(result[MISSING]>0);
    s->collision_by_hk[m.h-5][m.k-1]+=result[PAIRS];
    s->tested_by_hk[m.h-5][m.k-1]++;
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
static void matrix(FILE *out,uint64_t x[2][5]){
    fputc('[',out);
    for(Word h=0;h<2;h++){
        fputs(h?",[":"[",out);
        for(Word k=0;k<5;k++)fprintf(out,"%s%" PRIu64,k?",":"",x[h][k]);
        fputc(']',out);
    }
    fputc(']',out);
}
int main(int argc,char **argv){
    if(argc!=7)return 2;
    first=(Word)strtoul(argv[1],NULL,10);count=(Word)strtoul(argv[2],NULL,10);
    Word threads=(Word)strtoul(argv[3],NULL,10);double planned=strtod(argv[4],NULL);
    if(!count||!threads||threads>MAX_WORKERS)return 2;
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
        total.restricted_rank_losses+=per[i].restricted_rank_losses;
        total.collisions+=per[i].collisions;total.missed+=per[i].missed;
        for(Word j=0;j<VALID;j++)total.sum[j]+=per[i].sum[j];
        for(Word j=0;j<2;j++)for(Word k=0;k<5;k++){
            total.tested_by_hk[j][k]+=per[i].tested_by_hk[j][k];
            total.collision_by_hk[j][k]+=per[i].collision_by_hk[j][k];
        }
    }
    fprintf(out,"{\"first_id\":%u,\"count\":%u,\"threads\":%u,\"start_epoch\":%.6f,\"end_epoch\":%.6f,"
        "\"wall\":%.9f,\"cpu_seconds\":%.9f,\"positional_checksum64\":\"%016" PRIx64 "\",\"rank_losses\":%" PRIu64
        ",\"cases_with_collisions\":%" PRIu64 ",\"cases_with_missing_basis_bands\":%" PRIu64 ",\"sums\":[",
        first,count,threads,start,end,elapsed,cpu,total.checksum,total.restricted_rank_losses,total.collisions,total.missed);
    for(Word i=0;i<VALID;i++)fprintf(out,"%s%" PRIu64,i?",":"",total.sum[i]);
    fputs("],\"tested_by_hk\":",out);matrix(out,total.tested_by_hk);
    fputs(",\"collision_pairs_by_hk\":",out);matrix(out,total.collision_by_hk);
    fputs(",\"worker_cpu_seconds\":[",out);
    for(Word i=0;i<threads;i++)fprintf(out,"%s%.9f",i?",":"",per[i].worker_cpu_seconds);
    fputs("]}\n",out);
    return fclose(out)||fclose(telemetry);
}
