#include <inttypes.h>
#include <mach/host_info.h>
#include <mach/mach.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

enum { MAX=48, WORKERS=64, CATALOG=2*64*17 };
typedef struct {
    unsigned masks[MAX], weights[MAX], m, N, n, family, target;
} Cover;
typedef struct {
    uint64_t checksum,cases,overlap,tight_active,tight_rank;
    uint64_t groups[2][4];
    double work_seconds;
} Stats;
static uint32_t bank[CATALOG],first,count;
static atomic_uint next_block,active_threads;
static FILE *telemetry;

static uint32_t mix(uint32_t x) {
    x^=x>>16; x*=UINT32_C(0x7feb352d); x^=x>>15; x*=UINT32_C(0x846ca68b);
    return x^(x>>16);
}
static unsigned next(uint32_t *s){
    *s=*s*UINT32_C(1664525)+UINT32_C(1013904223);
    return *s>>12;
}
static void add(Cover*c,unsigned m,unsigned w){
    for(unsigned i=0;i<c->m;i++)if(c->masks[i]==m){c->weights[i]+=w;return;}
    if(c->m==MAX)abort();
    c->masks[c->m]=m;c->weights[c->m]=w;c->m++;
}
static void partition(Cover*c,unsigned id,unsigned w) {
    const uint32_t *p=bank+17*((c->n-3)*64+id);
    for(unsigned i=1;i<=p[0];i++)add(c,p[i],w);
}
static void cycles(Cover*c,uint32_t *s) {
    unsigned a[16],sizes[4]={3,5,3,5};
    for(unsigned i=0;i<c->N;i++)a[i]=i;
    for(unsigned i=c->N-1;i>0;i--){
        unsigned j=next(s)%(i+1),t=a[i];a[i]=a[j];a[j]=t;
    }
    unsigned start=0;
    for(unsigned j=0;j<c->N/4;j++){
        unsigned m=sizes[j];
        for(unsigned i=0;i<m;i++)add(c,(1u<<a[start+i])|(1u<<a[start+(i+1)%m]),1);
        start+=m;
    }
    if(start!=c->N)abort();
}
static Cover generate(unsigned id) {
    Cover c={0};
    c.n=3+((id>>2)&1u); c.family=id&3u;
    if(id<4)c.n=3+(id&1u);
    c.N=1u<<c.n;
    uint32_t seed=mix(id^UINT32_C(0xa2805901));
    if(id<2){
        for(unsigned j=0;j<c.n;j++){
            unsigned a=0;
            for(unsigned x=0;x<c.N;x++)if((x>>j)&1u)a|=1u<<x;
            add(&c,a,j+1);add(&c,((1u<<c.N)-1)^a,j+1);c.target+=j+1;
        }
    } else if(id<4 || c.family==2) {
        cycles(&c,&seed);c.target=2;
    } else if(c.family==3) {
        cycles(&c,&seed);partition(&c,next(&seed)%64,2);c.target=4;
    } else {
        for(unsigned j=0;j<2+c.family;j++){
            unsigned i=next(&seed)%64,w=next(&seed)%5+1;
            partition(&c,i,w);c.target+=w;
        }
    }
    return c;
}
static unsigned rank_gf2_add(unsigned b,unsigned pivot[4]){
    while(b){unsigned i=31u-__builtin_clz(b);
        if(!pivot[i]){pivot[i]=b;return 1;} b^=pivot[i];}
    return 0;
}
#define PRIME UINT32_C(2147483647)
static uint32_t reduce(uint64_t b){
    uint64_t a=(b&PRIME)+(b>>31);
    a=(a&PRIME)+(a>>31);
    return (uint32_t)(a>=PRIME ? a-PRIME:a);
}
static unsigned real_rank(const Cover*c) {
    uint32_t echelon[15][15]={{0}};
    unsigned rank=0;
    for(unsigned i=0;i<c->m;i++){
        uint32_t u[15];
        int base=c->masks[i]&1u;
        for(unsigned j=0;j<c->N-1;j++){
            int a=(int)((c->masks[i]>>(j+1))&1u)-base;
            u[j]=a<0 ? PRIME-1:(unsigned)a;
        }
        for(unsigned j=0;j<c->N-1;j++)if(u[j]){
            if(!echelon[j][j]){for(unsigned k=j;k<c->N-1;k++)echelon[j][k]=u[k];rank++;break;}
            uint32_t a=echelon[j][j],b=u[j];
            for(unsigned k=j;k<c->N-1;k++)
                u[k]=reduce((uint64_t)u[k]*a+(uint64_t)(PRIME-b)*echelon[j][k]);
        }
    }
    return rank;
}
static uint64_t checksum(unsigned id,unsigned code) {
    uint64_t x=((uint64_t)id<<24)|code;
    x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);
    return x^(x>>31);
}
static void audit(unsigned id,Stats*s) {
    Cover c=generate(id); unsigned largest=0,t=0,r=0,pivot[4]={0};
    for(unsigned x=0;x<c.N;x++){
        unsigned multiplicity=0,mass=0;
        for(unsigned i=0;i<c.m;i++)if((c.masks[i]>>x)&1u){multiplicity++;mass+=c.weights[i];}
        if(mass!=c.target)abort();
        if(multiplicity>largest)largest=multiplicity;
    }
    for(unsigned i=0;i<c.m;i++){
        unsigned m=c.masks[i],x0=__builtin_ctz(m);
        unsigned rank=0,p[4]={0};
        for(unsigned z=m;z;z&=z-1)rank+=rank_gf2_add(__builtin_ctz(z)^x0,p);
        if((1u<<rank)!=(unsigned)__builtin_popcount(m))abort();
        for(unsigned a=1;a<c.N;a++){
            unsigned yes=1;
            for(unsigned j=0;j<c.n;j++)if(__builtin_parity(a&p[j])){yes=0;break;}
            if(yes)t+=rank_gf2_add(a,pivot);
        }
    }
    r=real_rank(&c);
    if(t>c.m-largest || t>r)abort();
    unsigned code=c.m|(largest<<6)|(t<<12)|(r<<15)|((c.n-3)<<20)|(c.family<<21);
    s->checksum^=checksum(id,code);s->cases++;s->overlap+=largest>1;
    s->tight_active+=t==c.m-largest;s->tight_rank+=t==r;
    s->groups[c.n-3][c.family]++;
}
static double now(clockid_t id){struct timespec t;if(clock_gettime(id,&t))abort();return t.tv_sec+t.tv_nsec*1e-9;}
static void *work(void *arg) {
    Stats*s=arg; double a=now(CLOCK_THREAD_CPUTIME_ID);
    for(;;){
        unsigned i=atomic_fetch_add(&next_block,32);
        if(i>=count)break;
        unsigned hi=i+32; if(hi>count)hi=count;
        for(;i<hi;i++)audit(first+i,s);
    }
    s->work_seconds=now(CLOCK_THREAD_CPUTIME_ID)-a;
    atomic_fetch_sub(&active_threads,1);
    return NULL;
}
static host_cpu_load_info_data_t ticks(void){
    host_cpu_load_info_data_t t;mach_msg_type_number_t n=HOST_CPU_LOAD_INFO_COUNT;
    if(host_statistics(mach_host_self(),HOST_CPU_LOAD_INFO,(host_info_t)&t,&n)!=KERN_SUCCESS)abort();
    return t;
}
static void *monitor(void *unused){
    (void)unused;const struct timespec pause={0,120000000};
    double cpu=now(CLOCK_PROCESS_CPUTIME_ID),wall=now(CLOCK_MONOTONIC);
    host_cpu_load_info_data_t old=ticks();
    while(atomic_load(&active_threads)){
        nanosleep(&pause,NULL);
        double nc=now(CLOCK_PROCESS_CPUTIME_ID),nw=now(CLOCK_MONOTONIC);
        host_cpu_load_info_data_t fresh=ticks();
        unsigned long all=0,idle=fresh.cpu_ticks[CPU_STATE_IDLE]-old.cpu_ticks[CPU_STATE_IDLE];
        for(unsigned i=0;i<CPU_STATE_MAX;i++)all+=fresh.cpu_ticks[i]-old.cpu_ticks[i];
        fprintf(telemetry,"%.6f,%.9f,%.6f,%.6f,%u\n",now(CLOCK_REALTIME),nw-wall,
            all?100.0*(all-idle)/all:-1,100.0*(nc-cpu)/(nw-wall),atomic_load(&active_threads));
        fflush(telemetry);old=fresh;wall=nw;cpu=nc;
    }
    return NULL;
}
int main(int argc,char **argv){
    if(argc!=8)return 2;
    FILE*f=fopen(argv[1],"rb");if(!f||fread(bank,4,CATALOG,f)!=CATALOG)return 3;fclose(f);
    first=(unsigned)strtoul(argv[2],NULL,10);count=(unsigned)strtoul(argv[3],NULL,10);
    unsigned threads=(unsigned)strtoul(argv[4],NULL,10);
    double planned=strtod(argv[5],NULL);
    if(!count||!threads||threads>WORKERS)return 2;
    FILE*out=fopen(argv[6],"w");telemetry=fopen(argv[7],"w");if(!out||!telemetry)return 3;
    fputs("unix_epoch,sample_seconds,whole_machine_cpu_active_percent,a28_process_cpu_percent,active_workers\n",telemetry);
    while(planned>now(CLOCK_REALTIME)){const struct timespec p={0,5000000};nanosleep(&p,NULL);}
    atomic_init(&active_threads,threads);atomic_init(&next_block,0);
    pthread_t workers[WORKERS],sampler;Stats per[WORKERS]={0};
    double w=now(CLOCK_MONOTONIC),start=now(CLOCK_REALTIME),cpu=now(CLOCK_PROCESS_CPUTIME_ID);
    if(pthread_create(&sampler,NULL,monitor,NULL))abort();
    for(unsigned i=0;i<threads;i++)if(pthread_create(workers+i,NULL,work,per+i))abort();
    for(unsigned i=0;i<threads;i++)if(pthread_join(workers[i],NULL))abort();
    double end=now(CLOCK_REALTIME),elapsed=now(CLOCK_MONOTONIC)-w,cpusecs=now(CLOCK_PROCESS_CPUTIME_ID)-cpu;
    if(pthread_join(sampler,NULL))abort();
    Stats sum={0};
    for(unsigned i=0;i<threads;i++){
        sum.checksum^=per[i].checksum;sum.cases+=per[i].cases;sum.overlap+=per[i].overlap;
        sum.tight_active+=per[i].tight_active;sum.tight_rank+=per[i].tight_rank;
        for(unsigned j=0;j<2;j++)for(unsigned k=0;k<4;k++)sum.groups[j][k]+=per[i].groups[j][k];
    }
    if(sum.cases!=count)abort();
    fprintf(out,"{\"schema\":\"a28-r059-affine-exact-cover-cpu-v1\",\"first_id\":%u,\"count\":%u,"
      "\"threads\":%u,\"start_epoch\":%.6f,\"end_epoch\":%.6f,\"wall\":%.9f,\"cpu_seconds\":%.9f,"
      "\"cases\":%" PRIu64 ",\"overlapping_cases\":%" PRIu64 ",\"t_equals_m_minus_maximum\":%" PRIu64 ","
      "\"t_equals_centered_real_rank\":%" PRIu64 ",\"positional_checksum64\":\"%016" PRIx64 "\",\"by_dimension_and_family\":[",
      first,count,threads,start,end,elapsed,cpusecs,sum.cases,sum.overlap,sum.tight_active,sum.tight_rank,sum.checksum);
    for(unsigned i=0;i<2;i++)fprintf(out,"%s[%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 "]",
        i?",":"",sum.groups[i][0],sum.groups[i][1],sum.groups[i][2],sum.groups[i][3]);
    fputs("],\"worker_cpu_seconds\":[",out);
    for(unsigned i=0;i<threads;i++)fprintf(out,"%s%.9f",i?",":"",per[i].work_seconds);
    fputs("]}\n",out);
    return fclose(out)||fclose(telemetry);
}
