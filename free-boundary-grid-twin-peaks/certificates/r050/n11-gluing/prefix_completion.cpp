// Exact conditioned connected-growth completions in rectangles of size <=36.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using Word = std::uint64_t;
using Wide = __uint128_t;
constexpr Word ONE = 1;

struct U192 {
    std::array<Word,3> x{};
    explicit U192(Word value=0) : x{value,0,0} {}

    U192& operator+=(const U192& b) {
        Wide carry=0;
        for (int i=0;i<3;++i) {
            const Wide sum=Wide(x[i])+b.x[i]+carry;
            x[i]=Word(sum);
            carry=sum>>64;
        }
        if (carry) throw std::overflow_error("overwide sum");
        return *this;
    }
    U192& operator-=(const U192& b) {
        Wide borrow=0;
        for (int i=0;i<3;++i) {
            const Wide subtrahend=Wide(b.x[i])+borrow;
            const Wide minuend=x[i];
            x[i]=Word(minuend-subtrahend);
            borrow=minuend<subtrahend;
        }
        if (borrow) throw std::overflow_error("negative completion");
        return *this;
    }
    void times(Word factor) {
        Wide carry=0;
        for (Word& limb:x) {
            const Wide p=Wide(limb)*factor+carry;
            limb=Word(p);
            carry=p>>64;
        }
        if (carry) throw std::overflow_error("overwide product");
    }
    Word divide(Word d) {
        Word rem=0;
        constexpr Word mask=(ONE<<32)-1;
        for (int i=2;i>=0;--i) {
            const Word top=(rem<<32)|(x[i]>>32);
            const Word a=top/d;
            rem=top-a*d;
            const Word bottom=(rem<<32)|(x[i]&mask);
            const Word b=bottom/d;
            rem=bottom-b*d;
            x[i]=(a<<32)|b;
        }
        return rem;
    }
    std::string decimal() const {
        auto copy=*this;
        std::string out;
        do { out.push_back(char('0'+copy.divide(10))); }
        while (copy.x[0] || copy.x[1] || copy.x[2]);
        std::reverse(out.begin(),out.end());
        return out;
    }
};

class Grid {
public:
    Grid(int width,int height): w(width),h(height), n(w*h) {
        if (w<1 || w>6 || h<1 || h>6) throw std::invalid_argument("bad rectangle");
        full=(ONE<<n)-1;
        for (int v=0;v<n;++v) {
            const int x=v%w,y=v/w;
            if (x>0) adj[v]|=ONE<<(v-1);
            if (x+1<w) adj[v]|=ONE<<(v+1);
            if (y>0) adj[v]|=ONE<<(v-w);
            if (y+1<h) adj[v]|=ONE<<(v+w);
        }
    }

    Word neighbors(Word S) const {
        Word B=0;
        for (Word t=S;t;t&=t-1) B|=adj[std::countr_zero(t)];
        return B;
    }

    void validate(Word S) const {
        if (S==0 || (S&~full)) throw std::invalid_argument("prefix outside region or empty");
        Word reached=S&-S;
        for (;;) {
            const Word expanded=reached|(neighbors(reached)&S);
            if (expanded==S) return;
            if (expanded==reached) throw std::invalid_argument("prefix is disconnected");
            reached=expanded;
        }
    }

    int w,h,n;
    Word full=0;
    std::array<Word,36> adj{};
};

void enumerate(Word avail,Word S,const Grid& grid,std::vector<std::vector<Word>>& layers,Word& visited,Word limit) {
    if (++visited>limit) throw std::runtime_error("more states than planned");
    layers[std::popcount(S)].push_back(S);
    for (Word v=avail;v;v&=v-1) {
        const Word b=v&-v;
        enumerate((v^b)&~grid.adj[std::countr_zero(v)],S|b,grid,layers,visited,limit);
    }
}

struct Entry {
    Word key=0;
    U192 value{};
};
static_assert(sizeof(Entry)==32);

class Table {
public:
    explicit Table(Word slots,const U192& empty): data(slots), mask(slots-1), base(empty) {}
    const U192& lookup(Word key) const {
        if (!key) return base;
        for (Word ix=hash(key)&mask;;ix=(ix+1)&mask) {
            if (data[ix].key==key) return data[ix].value;
            if (!data[ix].key) throw std::logic_error("incomplete predecessor table");
        }
    }
    void insert(Word key,const U192& value) {
        Word ix=hash(key)&mask;
        while (data[ix].key) {
            if (data[ix].key==key) throw std::logic_error("duplicate independent set");
            ix=(ix+1)&mask;
        }
        data[ix]={key,value};
    }
private:
    static Word hash(Word k) {
        k=(k^(k>>30))*0xbf58476d1ce4e5b9ULL;
        k=(k^(k>>27))*0x94d049bb133111ebULL;
        return k^(k>>31);
    }
    std::vector<Entry> data;
    Word mask;
    U192 base;
};

struct Result {
    U192 count;
    Word states;
    Word predecessor_lookups;
    Word table_slots;
};

Result complete(const Grid& grid,Word S,Word planned) {
    grid.validate(S);
    const Word left=grid.full&~S;
    const Word forbidden=left&~grid.neighbors(S);
    const int n_left=std::popcount(left);
    U192 base(1);
    for (int j=2;j<=n_left;++j) base.times(j);
    std::vector<std::vector<Word>> layers(std::popcount(forbidden)+1);
    Word actual=0,queries=0;
    enumerate(forbidden,0,grid,layers,actual,planned);
    if (actual!=planned) throw std::logic_error("state count differs from independent plan");
    const Word slots=std::bit_ceil((4*actual+2)/3);
    Table table(slots,base);
    U192 positive=base,negative;
    for (std::size_t k=1;k<layers.size();++k) {
        for (Word I:layers[k]) {
            Word union_neighborhood=I;
            U192 value;
            for (Word q=I;q;q&=q-1) {
                value+=table.lookup(I^(q&-q));
                union_neighborhood|=grid.adj[std::countr_zero(q)]&left;
            }
            queries+=k;
            if (value.divide(std::popcount(union_neighborhood))) throw std::logic_error("nonintegral neighborhood recurrence");
            table.insert(I,value);
            (k%2==0 ? positive : negative)+=value;
        }
    }
    positive-=negative;
    return {positive,actual,queries,slots};
}

int main() try {
    int w,h;
    Word S,planned;
    while (std::cin>>w>>h>>S>>planned) {
        if (planned<1 || planned>1200000) throw std::invalid_argument("independent state budget refused");
        const Grid g(w,h);
        const Result r=complete(g,S,planned);
        std::cout << "{\"w\":"<<w<<",\"h\":"<<h<<",\"seen\":\""<<S<<"\",\"completion\":\""<<r.count.decimal()
                  <<"\",\"states\":"<<r.states<<",\"predecessor_lookups\":"<<r.predecessor_lookups
                  <<",\"table_slots\":"<<r.table_slots<<"}\n";
    }
    if (!std::cin.eof()) throw std::invalid_argument("malformed input");
    return 0;
} catch (const std::exception& e) {
    std::cerr<<e.what()<<'\n';
    return 2;
}
