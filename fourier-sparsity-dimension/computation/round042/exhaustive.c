// Complete balanced depth-three trees on six specified two-variable Boolean leaves.
#include <stdio.h>
#include <stdlib.h>
static const unsigned leaves[]={0,15,3,5,1,2};
static void fwht(int *a,int n){
    for(int w=1;w<n;w*=2)for(int p=0;p<n;p+=2*w)for(int j=p;j<p+w;j++){
        int x=a[j],y=a[j+w];a[j]=x+y;a[j+w]=x-y;
    }
}
static int rank(const int *in,int n,int bits){
    int basis[10]={0},r=0;
    for(int i=0;i<n;i++)if(in[i]){
        int v=i;
        for(int b=bits-1;b>=0;b--)if(v&(1<<b)){
            if(!basis[b]){basis[b]=v;r++;break;}
            v^=basis[b];
        }
    }
    return r;
}
int main(int argc,char**argv){
    if(argc!=3)return 2;
    int first=atoi(argv[1]),step=atoi(argv[2]);
    if(step<1||first<0||first>=step)return 2;
    int route[128];
    for(int z=0;z<128;z++){
        int v=1;
        while(v<8)v=2*v+((z>>(v-1))&1);
        route[z]=v;
    }
    int count=0,tight=0;
    for(int id=first;id<1679616;id+=step){
        int sum[16][4]={0},code=id,table[512],select[128]={0},tags[4]={0};
        int predicted=0,mark=0;
        for(int v=8;v<16;v++,code/=6)for(int y=0;y<4;y++)
            sum[v][y]=1-2*((leaves[code%6]>>y)&1);
        for(int v=7;v>=1;v--){
            int diff[4];
            for(int y=0;y<4;y++){
                sum[v][y]=sum[2*v][y]+sum[2*v+1][y];
                diff[y]=sum[2*v][y]-sum[2*v+1][y];
            }
            fwht(diff,4);
            int terms=0;
            for(int y=0;y<4;y++)terms+=diff[y]!=0;
            predicted+=(v>=4?4:v>=2?2:1)*terms;
            if(terms)for(int a=v;a;a/=2)mark|=1<<(a-1);
        }
        fwht(sum[1],4);
        for(int y=0;y<4;y++)predicted+=sum[1][y]!=0;
        for(int z=0;z<128;z++)for(int y=0;y<4;y++)
            table[z*4+y]=sum[route[z]][y];
        fwht(table,512);
        int s=0;
        for(int i=0;i<512;i++)if(table[i]){s++;select[i/4]=1;tags[i%4]=1;}
        int n=rank(select,128,7),d=rank(table,512,9);
        if(s!=predicted||n!=__builtin_popcount((unsigned)mark)||(n+1)*(n+1)>4*s||
            d>n+rank(tags,4,2))return 1;
        tight+=(n+1)*(n+1)==4*s;count++;
    }
    printf("%d %d %d\n",first,count,tight);
}
