// Enumerate set partitions of F_2^4 into at most four nonempty cells.
#include <array>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>
using Mask = unsigned;
constexpr Mask full = 65535;
struct Cell { Mask spectrum; int size; bool flat; };
std::array<Cell,65536> cells;
std::array<int,65536> ranks;
struct Counts { uint64_t n=0, sharp=0, nonflat=0, failures=0; };
void test(const std::array<Mask,4>& blocks,int b, Counts& out) {
    int total=0, badCost=0, badCount=0;
    Mask all=0, bad=0;
    for (int i=0;i<b;++i) {
        auto c=cells[blocks[i]];
        all|=c.spectrum; total+=c.size;
        if (!c.flat) {bad|=c.spectrum; badCost+=c.size; ++badCount;}
    }
    int l=std::bit_floor(unsigned(b)), delta=total-l*(3*b-2*l);
    int g=b-badCount;
    ++out.n; out.sharp+=(delta==0); out.nonflat+=(badCount!=0);
    if(delta<0 || badCost>2*delta || ranks[all]>ranks[bad]+(g?g-1:0)
       || ranks[all]>b-1+2*delta) ++out.failures;
}
void visit(Mask rest,std::array<Mask,4>& blocks,int used,Counts& out) {
    blocks[used]=rest; test(blocks,used+1,out);
    if(used==3) return;
    Mask free=rest&(rest-1), first=rest^free;
    for(Mask t=free;t;t=(t-1)&free) {
        blocks[used]=first|(free^t); visit(t,blocks,used+1,out);
    }
}
int main(int argc,char** argv) {
    if(argc!=4) return 2;
    int a=std::atoi(argv[1]), b=std::atoi(argv[2]), threads=std::atoi(argv[3]);
    if(a<0 || a>=b || b>32768 || threads<1) return 2;
    for(Mask mask=1;mask<=full;++mask) {
        Mask spec=0;
        for(int u=0;u<16;++u) {
            int sum=0;
            for(int x=0;x<16;++x) if((mask>>x)&1) sum+=1-2*(std::popcount(unsigned(u&x))%2);
            if(sum) spec|=1u<<u;
        }
        int size=std::popcount(spec);
        cells[mask]={spec,size,size*std::popcount(mask)==16};
        Mask span=1;
        for(int u=1;u<16;++u) if(((mask>>u)&1) && !((span>>u)&1)) {
            Mask old=span;
            for(int x=0;x<16;++x) if((old>>x)&1) span|=1u<<(x^u);
            ++ranks[mask];
        }
    }
    std::vector<Counts> counts(threads);
    std::vector<std::thread> jobs;
    for(int w=0;w<threads;++w) jobs.emplace_back([&,w] {
        std::array<Mask,4> blocks{};
        for(int j=a+w;j<b;j+=threads) {
            blocks[0]=1u|(unsigned(j)<<1);
            Mask rest=full^blocks[0];
            if(rest) visit(rest,blocks,1,counts[w]); else test(blocks,1,counts[w]);
        }
    });
    Counts result;
    for(auto& t:jobs)t.join();
    for(auto c:counts){result.n+=c.n;result.sharp+=c.sharp;result.nonflat+=c.nonflat;result.failures+=c.failures;}
    std::printf("{\"first\":%d,\"end\":%d,\"threads\":%d,\"partitions\":%llu,\"minimum_cost\":%llu,\"with_nonflat\":%llu,\"failures\":%llu}\n",a,b,threads,(unsigned long long)result.n,(unsigned long long)result.sharp,(unsigned long long)result.nonflat,(unsigned long long)result.failures);
    return result.failures!=0;
}
