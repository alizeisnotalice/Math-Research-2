// Count fundamental cuts for every labeled spanning tree of a triangle necklace.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
static int analyze(int m,unsigned code,uint64_t *two,uint64_t *four){
    const unsigned rim=(1u<<m)-1;
    const unsigned broken=code%(unsigned)m;
    code/=m;
    unsigned choice=code%2;
    code/=2;
    unsigned adj[30]={0};
    for(int i=0;i<m;i++){
        int j=(i+1)%m;
        unsigned edges;
        if((unsigned)i==broken)edges=choice?2:4;
        else {edges=7u^(1u<<(code%3));code/=3;}
        int a[3]={i,i,j},b[3]={j,m+i,m+i};
        for(int k=0;k<3;k++)if(edges&(1u<<k)){
            adj[a[k]]|=1u<<b[k];adj[b[k]]|=1u<<a[k];
        }
    }
    int order[30]={0},parent[30]={0},n=1;
    unsigned seen=1,sub[30]={0};
    for(int k=0;k<n;k++){
        int v=order[k];
        unsigned children=adj[v]&~seen;
        seen|=children;
        while(children){
            int u=__builtin_ctz(children);
            children&=children-1;parent[u]=v;order[n++]=u;
        }
    }
    if(n!=2*m)return 0;
    for(int k=n-1;k>0;k--){
        int v=order[k];sub[v]|=1u<<v;sub[parent[v]]|=sub[v];
        unsigned a=sub[v]&rim,b=(a>>1)|((a&1u)<<(m-1)),c=sub[v]>>m;
        int cut=__builtin_popcount(a^b)+__builtin_popcount(a^c)+__builtin_popcount(b^c);
        if(cut==2)++*two;else if(cut==4)++*four;else return 0;
    }
    return 1;
}
int main(int argc,char **argv){
    if(argc!=4)return 2;
    int m=atoi(argv[1]),shard=atoi(argv[2]),step=atoi(argv[3]);
    if(m<3||m>14||shard<0||step<1||shard>=step)return 2;
    unsigned total=2*m;
    for(int i=1;i<m;i++)total*=3;
    uint64_t two=0,four=0;unsigned count=0;
    for(unsigned i=shard;i<total;i+=step){
        if(!analyze(m,i,&two,&four))return 1;
        count++;
    }
    printf("%d %u %llu %llu\n",shard,count,(unsigned long long)two,(unsigned long long)four);
}
