// Exact Boolean meet in the middle for negative reciprocal C modulo J^2.
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int p=163;
using Sig=std::array<int,14>;
using Key=__uint128_t;
constexpr Key sentinel=~Key{0};

int residue(int value) {
    int r=value%p;
    return r<0 ? r+p : r;
}

void add(Sig& left, const Sig& right, bool negative) {
    for (unsigned j=0;j<14;++j) {
        const int x=left[j]+(negative ? -right[j] : right[j]);
        left[j]=x+(x<0)*p-(x>=p)*p;
    }
}

Key pack(const Sig& v) {
    Key key=0;
    for (int r:v) key=key*p+static_cast<uint64_t>(r);
    return key;
}

uint64_t mix(Key packed) {
    uint64_t key=static_cast<uint64_t>(packed)^std::rotl(static_cast<uint64_t>(packed>>64),23);
    key=(key^(key>>30))*0xbf58476d1ce4e5b9ULL;
    key=(key^(key>>27))*0x94d049bb133111ebULL;
    return key^(key>>31);
}

std::vector<Sig> columns(int degree,int q) {
    std::array<int,8> j{q,117,115,68,106,121,6,1};
    std::array<int,15> modulus{};
    for (unsigned i=0;i<8;++i)
        for (unsigned k=0;k<8;++k)
            modulus[i+k]=residue(modulus[i+k]+j[i]*j[k]);
    std::vector<Sig> basis(degree+1), result(degree+1);
    basis[0][0]=1;
    for (int k=1;k<=degree;++k) {
        int overflow=basis[k-1][13];
        for (int j=0;j<14;++j)
            basis[k][j]=residue((j==0 ? 0 : basis[k-1][j-1])-overflow*modulus[j]);
    }
    for (int i=0;i<=degree;++i) {
        for (int j=0;j<14;++j)
            result[i][j]=residue((i%2==0 ? 1 : -1)*basis[degree-i][j]);
    }
    return result;
}

class Table {
public:
    explicit Table(uint64_t size): keys(size,sentinel), masks(size), mask(size-1) {}

    void insert(Key key,uint32_t bits) {
        uint64_t slot=mix(key)&mask;
        while (keys[slot]!=sentinel) slot=(slot+1)&mask;
        keys[slot]=key;
        masks[slot]=bits;
    }

    void emit_matches(Key key,uint64_t prefix,bool& first) const {
        uint64_t slot=mix(key)&mask;
        while (keys[slot]!=sentinel) {
            if (keys[slot]==key) {
                if (!first) std::cout << ',';
                first=false;
                std::cout << (prefix|masks[slot]);
            }
            slot=(slot+1)&mask;
        }
    }
private:
    std::vector<Key> keys;
    std::vector<uint32_t> masks;
    uint64_t mask;
};
}

int main(int argc,char** argv) {
    if (argc!=3 && argc!=4)
        throw std::invalid_argument("usage: degree q [fixture_mask]");
    const int degree=std::stoi(argv[1]);
    const int q=std::stoi(argv[2]);
    const bool fixture=argc==4;
    if (degree<13 || degree>64 || q<0 || q>=p || (fixture && degree>30))
        throw std::invalid_argument("degree, q or fixture range");
    const unsigned total=degree-12;
    const unsigned left_bits=total/2;
    const uint64_t left_size=uint64_t{1}<<left_bits;
    const uint64_t right_size=uint64_t{1}<<(total-left_bits);
    const auto col=columns(degree,q);

    Sig target{};
    if (fixture) {
        uint64_t known=std::stoull(argv[3]);
        if (known>=(uint64_t{1}<<total))
            throw std::invalid_argument("fixture mask out of range");
        for (unsigned i=0;i<total;++i)
            if ((known>>i)&1) add(target,col[12+i],false);
    } else {
        for (unsigned i=0;i<12;++i)
            if (i!=7) add(target,col[i],true);
        add(target,col[degree],true);
    }

    Table table(2*left_size);
    Sig left{};
    uint64_t previous=0;
    for (uint64_t i=0;i<left_size;++i) {
        uint64_t gray=i^(i>>1);
        if (i!=0) {
            uint64_t change=gray^previous;
            add(left,col[12+std::countr_zero(change)],(gray&change)==0);
        }
        table.insert(pack(left),static_cast<uint32_t>(gray));
        previous=gray;
    }

    std::cout << "{\"degree\":" << degree << ",\"q\":" << q
              << ",\"mode\":\"" << (fixture ? "fixture" : "boundary") << "\",\"left_patterns\":"
              << left_size << ",\"right_patterns\":" << right_size << ",\"matches\":[";
    previous=0;
    bool first=true;
    for (uint64_t i=0;i<right_size;++i) {
        uint64_t gray=i^(i>>1);
        if (i!=0) {
            uint64_t change=gray^previous;
            add(target,col[12+left_bits+std::countr_zero(change)],(gray&change)!=0);
        }
        table.emit_matches(pack(target),gray<<left_bits,first);
        previous=gray;
    }
    std::cout << "]}\n";
}
