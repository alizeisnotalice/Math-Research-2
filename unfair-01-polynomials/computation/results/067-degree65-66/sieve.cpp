// Exact coefficient-slice Euclidean sieve in the ordinary x variable.
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
constexpr int maximum_degree=66;
using Poly=std::array<int,maximum_degree+1>;

constexpr int residue(int x) {
    const int r=x%p;
    return r<0 ? r+p : r;
}

constexpr auto inverse_table() {
    std::array<int,p> a{};
    a[1]=1;
    for (int i=2;i<p;++i) a[i]=p-(p/i)*a[p%i]%p;
    return a;
}

int quotient_strip(int* left,int n,const int* right,int m) {
    constexpr auto inv=inverse_table();
    const int d=n-m;
    const int scale=inv[right[0]];
    // In the usual degree-drop-one case, solve both quotient coefficients first.
    // There is then only one modular reduction per residual coefficient.
    if (d==1) {
        const int u=left[0]*scale%p;
        const int v=residue(left[1]-u*right[1])*scale%p;
        for (int i=2;i<=m;++i)
            left[i]=(left[i]+2*p*p-u*right[i]-v*right[i-1])%p;
        left[n]=(left[n]+p*p-v*right[m])%p;
    } else {
        Poly quotient{};
        for (int i=0;i<=d;++i) {
            int top=left[i];
            for (int j=std::max(0,i-m);j<i;++j) top-=quotient[j]*right[i-j];
            quotient[i]=residue(top)*scale%p;
        }
        for (int i=d+1;i<=n;++i) {
            int value=left[i];
            for (int j=std::max(0,i-m);j<=d;++j) value-=quotient[j]*right[i-j];
            left[i]=residue(value);
        }
    }
    int start=d+1;
    while (start<=n && left[start]==0) ++start;
    return start;
}

int gcd_degree(Poly a,Poly b,int n,int threshold) {
    int* left=a.data();
    int* right=b.data()+1;
    int m=n-1;
    while (m>=0 && *right==0) {
        ++right;
        --m;
    }
    while (m>=0) {
        if (!m || m<threshold) return m;
        const int stripped=quotient_strip(left,n,right,m);
        n-=stripped;
        left+=stripped;
        std::swap(left,right);
        std::swap(n,m);
    }
    return n;
}

std::vector<int> exponents(const std::string& str,int n) {
    std::vector<int> positions;
    std::size_t at=0;
    if (str=="-") return positions;
    while (at<str.size()) {
        const auto end=str.find(',',at);
        const auto tok=str.substr(at,end-at);
        std::size_t read=0;
        const int i=std::stoi(tok,&read);
        if (i<0 || i>n || read!=tok.size() || (!positions.empty() && i<=positions.back()))
            throw std::invalid_argument("expected distinct increasing exponents within degree");
        positions.push_back(i);
        if (end==std::string::npos) break;
        at=end+1;
    }
    return positions;
}

void set_bit(Poly& C,Poly& D,int n,int k,int value) {
    C[n-k]=value;
    // Logarithmic derivative x*C' minus n*C.  gcd(C,x)=1 since C(0)=1.
    D[n-k]=value ? residue(k-n) : 0;
}

void fixed_rows(Poly& C,Poly& D,const std::vector<int>& ones,int n) {
    if (ones.empty() || ones.front()!=0 || ones.back()!=n)
        throw std::invalid_argument("fixed terms must contain both endpoints");
    for (const int k:ones) set_bit(C,D,n,k,1);
}
}

int main(int argc,char** argv) {
    if (argc==4 && std::string(argv[1])=="--inspect") {
        const int n=std::stoi(argv[2]);
        if (n<8 || n>maximum_degree) throw std::invalid_argument("degree");
        Poly a{},b{};
        fixed_rows(a,b,exponents(argv[3],n),n);
        std::cout << gcd_degree(a,b,n,0) << '\n';
        return 0;
    }
    if (argc!=6) throw std::invalid_argument("degree one_csv free_csv begin end");
    const int n=std::stoi(argv[1]);
    if (n<8 || n>maximum_degree) throw std::invalid_argument("degree");
    const auto ones=exponents(argv[2],n);
    const auto free=exponents(argv[3],n);
    if (free.size()>26) throw std::invalid_argument("more than 26 free bits");
    for (const int i:free)
        if (std::binary_search(ones.begin(),ones.end(),i))
            throw std::invalid_argument("fixed/free overlap");
    const auto begin=std::stoull(argv[4]),end=std::stoull(argv[5]);
    if (begin>=end || end>(uint64_t{1}<<free.size()))
        throw std::invalid_argument("range");
    Poly a{},b{};
    fixed_rows(a,b,ones,n);
    uint64_t gray=begin^(begin>>1);
    for (unsigned i=0;i<free.size();++i)
        if ((gray>>i)&1) set_bit(a,b,n,free[i],1);
    std::cout << "{\"degree\":" << n << ",\"begin\":" << begin << ",\"end\":" << end << ",\"hits\":[";
    bool first=true;
    for (uint64_t i=begin;i<end;++i) {
        const auto next=i^(i>>1);
        if (next!=gray) {
            const unsigned bit=std::countr_zero(next^gray);
            set_bit(a,b,n,free[bit],(next>>bit)&1);
        }
        gray=next;
        const int d=gcd_degree(a,b,n,7);
        if (d<7) continue;
        if (!first) std::cout << ',';
        first=false;
        std::cout << '[' << i << ',' << gray << ',' << d << ']';
    }
    std::cout << "],\"count\":" << end-begin << "}\n";
}
