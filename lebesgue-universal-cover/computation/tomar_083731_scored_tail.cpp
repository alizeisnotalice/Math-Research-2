struct ScoredTraversal {
    Tester& fan;
    const std::uint8_t* codes;
    const std::vector<int>& stop;
    const std::vector<int>& prefix;
    const L* radii;
    const L* radial_floors;
    std::uint8_t* pending;
    L* raw_scores;
    L* local_scores;
    L* sums;
    int max_depth, max_nodes;

    bool probe(Box& lo, Box& hi, int depth, L& nodes, L& minimum) {
        if(nodes==max_nodes)return false;
        ++nodes;
        sums[1]=std::max(sums[1],L(depth));
        if(fan.check(lo,hi)==1){
            minimum=std::min(minimum,fan.current_floor);
            ++sums[4];
            return true;
        }
        if(depth==max_depth)return false;
        const int k=fan.split_axis;
        const L l=lo[k],h=hi[k];
        if(h-l<2)return false;
        const L m=l+(h-l)/2;
        hi[k]=m;
        const bool left=probe(lo,hi,depth+1,nodes,minimum);
        hi[k]=h;
        if(!left)return false;
        lo[k]=m;
        const bool right=probe(lo,hi,depth+1,nodes,minimum);
        lo[k]=l;
        return right;
    }

    L exterior_cap(const Box& lo,const Box& hi,const Box& clipped_lo,const Box& clipped_hi) const{
        L ans=std::numeric_limits<L>::max();
        for(int b=0;b<4;++b){
            const int j=3*b;
            for(int k=j;k<=j+1;++k)
                if(lo[k]!=clipped_lo[k] || hi[k]!=clipped_hi[k])
                    ans=std::min(ans,radial_floors[b]);
        }
        return ans;
    }

    void emit(int at,int count,L inside,L exterior){
        const L full=inside?std::min(inside,exterior):exterior;
        if(full<fan.goal)std::abort();
        for(int j=at;j<stop[at];++j){
            if(prefix[j+1]>prefix[j]){
                raw_scores[j]=full;
                local_scores[j]=inside;
            }
        }
        sums[2]+=count;
    }

    bool walk(int at,Box& lo,Box& hi,int depth){
        if(depth>256)return false;
        const int count=prefix[stop[at]]-prefix[at];
        if(!count)return true;
        const int code=codes[at];
        if(code>=253 || (count>=64 && depth%5==0)){
            Box l=lo,h=hi;
            if(!clip_domain(l,h,radii)){
                emit(at,count,0,*std::min_element(radial_floors,radial_floors+4));
                return true;
            }
            if(code>=253){
                L minimum=std::numeric_limits<L>::max();
                L visits=0;
                const bool ok=probe(l,h,0,visits,minimum);
                sums[0]+=visits;
                if(ok)emit(at,1,minimum,exterior_cap(lo,hi,l,h));
                else{
                    pending[at]=1;
                    ++sums[3];
                    sums[5]+=(visits==max_nodes);
                }
                return true;
            }
            ++sums[0];
            if(fan.check(l,h)==1){
                emit(at,count,fan.current_floor,exterior_cap(lo,hi,l,h));
                return true;
            }
        }

        const L a=lo[code], b=hi[code];
        if((a+b)%2)return false;
        const L m=(a+b)/2;
        hi[code]=m;
        if(!walk(at+1,lo,hi,depth+1))return false;
        hi[code]=b;lo[code]=m;
        if(!walk(stop[at+1],lo,hi,depth+1))return false;
        lo[code]=a;
        return true;
    }
};

extern "C" int adaptive_scored_trim(
    const std::uint8_t* codes,int n,const L* lower,const L* upper,const L* orig_hi,
    const std::uint8_t* active,std::uint8_t* pending,L* stats,const std::int32_t* table,
    const std::int32_t* sizes,L goal,int depth,int calls,const L* radii,
    const L* radial_floors,L* raw_scores,L* local_scores){
    if(n<1 || n>(1<<29) || depth<0 || depth>64 || calls<1 || calls>(1<<24))return 1;
    for(int j=0;j<11;++j)
        if(lower[j]<=-Fine || upper[j]>=Fine || lower[j]>=upper[j])return 2;
    for(int j=0;j<4;++j)
        if(radii[j]<=0 || radial_floors[j]<goal)return 3;
    std::vector<int> prefix(n+1),stop(n);
    for(int j=0;j<n;++j){
        if(active[j]>1 || (active[j] && codes[j]<253))return 4;
        prefix[j+1]=prefix[j]+active[j];
    }
    for(int j=n-1;j>=0;--j){
        if(codes[j]>=253)stop[j]=j+1;
        else{
            if(codes[j]>=11 || j+1>=n || stop[j+1]>=n)return 5;
            stop[j]=stop[stop[j+1]];
        }
    }
    if(stop[0]!=n)return 6;
    Tester fan(table,sizes,orig_hi,goal);
    ScoredTraversal r{fan,codes,stop,prefix,radii,radial_floors,pending,raw_scores,local_scores,
                      stats,depth,calls};
    Box lo,hi;
    std::copy(lower,lower+11,lo.begin()); std::copy(upper,upper+11,hi.begin());
    return r.walk(0,lo,hi,0)?0:7;
}
