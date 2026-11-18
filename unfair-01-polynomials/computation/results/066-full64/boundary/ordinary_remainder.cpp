// Exact projected join in the ordinary x coordinate; every hit is filtered by all 14 residues.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int p=163;
constexpr int width=14;
constexpr int projection_width=5;
constexpr std::array<int,projection_width> axes{0,3,6,9,12};
constexpr uint64_t code_count=uint64_t{p}*p*p*p*p;
constexpr unsigned bucket_shift=18;
constexpr uint64_t buckets=(code_count+(uint64_t{1}<<bucket_shift)-1)>>bucket_shift;
using Full=std::array<int,width>;
using Short=std::array<int,projection_width>;
using Modulus=std::array<int,width+1>;

int reduced(int n) {
    const int t=n%p;
    return t<0 ? t+p : t;
}

int inv(int a) {
    for (int b=1;b<p;++b)
        if ((a*b)%p==1) return b;
    throw std::invalid_argument("nonunit q");
}

Modulus ordinary_modulus(int q) {
    // K(x)=x^7 J(-1/x)/q, with coefficients in ascending order.
    const std::array<int,8> j{q,117,115,68,106,121,6,1};
    std::array<int,8> k{};
    const int normal=inv(q);
    for (unsigned i=0;i<8;++i)
        k[7-i]=reduced((i%2 ? -1 : 1)*normal*j[i]);
    Modulus out{};
    for (unsigned i=0;i<8;++i)
        for (unsigned h=0;h<8;++h)
            out[i+h]=reduced(out[i+h]+k[i]*k[h]);
    if (out.back()!=1) throw std::logic_error("nonmonic ordinary modulus");
    return out;
}

Full remainder(int exponent,const Modulus& modulus) {
    std::array<int,65> row{};
    row[exponent]=1;
    for (int i=exponent;i>=width;--i) {
        const int c=row[i];
        if (c==0) continue;
        for (unsigned j=0;j<width;++j)
            row[i-width+j]=reduced(row[i-width+j]-c*modulus[j]);
        row[i]=0;
    }
    Full out{};
    std::copy_n(row.begin(),width,out.begin());
    return out;
}

Short project(const Full& row) {
    Short out{};
    for (unsigned i=0;i<axes.size();++i) out[i]=row[axes[i]];
    return out;
}

Short shifted(const Short& a,const Short& b,bool subtract=false) {
    Short out{};
    for (unsigned i=0;i<out.size();++i) {
        const int x=a[i]+(subtract ? -b[i] : b[i]);
        out[i]=x+(x<0)*p-(x>=p)*p;
    }
    return out;
}

uint64_t key(const Short& row) {
    uint64_t out=0;
    for (const int x:row) out=out*p+x;
    return out;
}

std::vector<Short> subsets(const std::vector<Short>& input,unsigned start,unsigned count) {
    std::vector<Short> sums(uint64_t{1}<<count);
    for (unsigned mask=1;mask<sums.size();++mask) {
        const unsigned low=std::countr_zero(mask);
        sums[mask]=shifted(sums[mask&(mask-1)],input[start+low]);
    }
    return sums;
}

template<class Visit>
void binary_half(const std::vector<Short>& cols, Visit visit) {
    const unsigned lower=cols.size()/2;
    const auto a=subsets(cols,0,lower);
    const auto b=subsets(cols,lower,cols.size()-lower);
    for (unsigned upper=0;upper<b.size();++upper)
        for (unsigned low=0;low<a.size();++low)
            visit((uint64_t{upper}<<lower)|low,shifted(a[low],b[upper]));
}

bool confirmed(const Full& wanted,const std::vector<Full>& cols,uint64_t mask) {
    Full difference=wanted;
    while (mask) {
        const unsigned i=std::countr_zero(mask);
        mask&=mask-1;
        for (unsigned j=0;j<difference.size();++j)
            difference[j]-=cols[i][j];
    }
    return std::all_of(difference.begin(),difference.end(),[](int x){return x%p==0;});
}

template<class R>
void print_array(const R& row) {
    std::cout << '[';
    bool first=true;
    for (const auto x:row) {
        if (!first) std::cout << ',';
        first=false;
        std::cout << x;
    }
    std::cout << ']';
}
}

int main(int argc,char** argv) {
    if (argc<3 || argc>4) throw std::invalid_argument("usage: degree q [--tight | --fixture=mask]");
    const int degree=std::stoi(argv[1]),q=std::stoi(argv[2]);
    const std::string mode=argc==4 ? argv[3] : "";
    const bool tight=mode=="--tight", fixture=mode.starts_with("--fixture=");
    if (degree<13 || degree>64 || q<1 || q>=p || (tight && degree!=64) ||
        (fixture && degree>30) || (!mode.empty() && !fixture && !tight))
        throw std::invalid_argument("degree, q or mode");
    std::array<int,65> bits{};
    bits.fill(-1);
    if (!fixture) {
        for (int j=0;j<12;++j) bits[j]=(j!=7);
        bits[degree]=1;
        if (tight) {
            for (int j:{31,32,33,34,35,38,58}) bits[j]=1;
            for (int j:{60,61,62,63}) bits[j]=0;
        }
    }

    const auto mod=ordinary_modulus(q);
    std::vector<Full> all(degree+1);
    for (int k=0;k<=degree;++k) all[k]=remainder(k,mod);
    std::vector<int> positions;
    std::vector<Full> variable;
    Full wanted{};
    uint64_t fixed_mask=0;
    const uint64_t fixture_mask=fixture ? std::stoull(mode.substr(10)) : 0;
    if (fixture && fixture_mask>=(uint64_t{1}<<(degree-12))) throw std::invalid_argument("fixture mask");
    for (int k=0;k<=degree;++k) {
        if (k>=12 && k<degree && bits[k]<0) {
            positions.push_back(k);
            variable.push_back(all[k]);
            if (fixture && ((fixture_mask>>(k-12))&1))
                for (unsigned j=0;j<width;++j) wanted[j]=reduced(wanted[j]+all[k][j]);
        } else if (!fixture && bits[k]==1) {
            for (unsigned j=0;j<width;++j) wanted[j]=reduced(wanted[j]-all[k][j]);
            if (k>=12 && k<degree) fixed_mask|=uint64_t{1}<<(k-12);
        }
    }
    const unsigned stored=positions.size()/2;
    const unsigned scanned=positions.size()-stored;
    const uint64_t low_mask=(uint64_t{1}<<stored)-1;
    const uint64_t stored_count=uint64_t{1}<<stored;
    const uint64_t scan_count=uint64_t{1}<<scanned;
    std::vector<Short> scan,store;
    for (unsigned i=0;i<variable.size();++i)
        (i<scanned ? scan : store).push_back(project(variable[i]));

    // Materialize the high exponents in order-independent sorted buckets.
    std::vector<uint64_t> table(stored_count);
    std::vector<uint32_t> limits(buckets+1);
    binary_half(store,[&](uint64_t, const Short& sum){++limits[(key(sum)>>bucket_shift)+1];});
    for (uint64_t i=0;i<buckets;++i) limits[i+1]+=limits[i];
    if (limits.back()!=stored_count) throw std::logic_error("missing stored masks");
    auto cursor=limits;
    binary_half(store,[&](uint64_t mask,const Short& sum){
        const auto k=key(sum);
        table[cursor[k>>bucket_shift]++]=(k<<stored)|mask;
    });
    for (uint64_t i=0;i<buckets;++i) {
        if (cursor[i]!=limits[i+1]) throw std::logic_error("incomplete bucket");
        if (limits[i+1]-limits[i]>1)
            std::sort(table.begin()+limits[i],table.begin()+limits[i+1]);
    }
    std::vector<uint64_t> matches,projected_masks;
    uint64_t projected_pairs=0,scanned_count=0;
    const auto target=project(wanted);
    binary_half(scan,[&](uint64_t mask,const Short& sum){
        ++scanned_count;
        const auto k=key(shifted(target,sum,true));
        const auto bucket=k>>bucket_shift;
        const auto end=table.cbegin()+limits[bucket+1];
        auto start=std::lower_bound(table.cbegin()+limits[bucket],end,k<<stored);
        for (;start!=end && (*start>>stored)==k;++start) {
            ++projected_pairs;
            const auto candidate=mask|((*start&low_mask)<<scanned);
            uint64_t old_index=fixed_mask;
            for (unsigned i=0;i<positions.size();++i)
                if ((candidate>>i)&1) old_index|=uint64_t{1}<<(positions[i]-12);
            projected_masks.push_back(old_index);
            if (confirmed(wanted,variable,candidate)) matches.push_back(old_index);
        }
    });
    if (scanned_count!=scan_count) throw std::logic_error("missing scanned masks");
    std::sort(matches.begin(),matches.end());
    std::sort(projected_masks.begin(),projected_masks.end());
    if (projected_masks.size()!=projected_pairs ||
         std::adjacent_find(projected_masks.begin(),projected_masks.end())!=projected_masks.end())
        throw std::logic_error("incomplete or duplicate projection candidates");
    if (std::adjacent_find(matches.begin(),matches.end())!=matches.end())
        throw std::logic_error("duplicate full mask");
    std::cout << "{\"degree\":" << degree << ",\"q\":" << q << ",\"mode\":\""
              << (tight ? "tight" : fixture ? "fixture" : "same_range")
              << "\",\"stored_patterns\":" << stored_count << ",\"scanned_patterns\":" << scan_count
              << ",\"projected_pairs\":" << projected_pairs << ",\"fixed_mask\":" << fixed_mask
              << ",\"positions\":";
    print_array(positions);
    std::cout << ",\"modulus\":";
    print_array(mod);
    std::cout << ",\"target_remainder\":";
    print_array(wanted);
    std::cout << ",\"column_remainders\":[";
    for (unsigned i=0;i<all.size();++i) {
        if (i) std::cout << ',';
        print_array(all[i]);
    }
    std::cout << "],\"projected_masks\":";
    print_array(projected_masks);
    std::cout << ",\"matches\":";
    print_array(matches);
    std::cout << "}\n";
}
