#define main production_screen_main
#include "s3_tail_screen.cpp"
#undef main

int main() {
    // 2 P where P=z^7(z-r), r=+-2, at a strict boundary. Single
    // nonzero lower high coefficient; generic kernel synthetic fixture.
    Group g{};
    g.D = 5;
    const auto enclosed = [](const Poly& p) {
        std::array<Interval,9> a{};
        for (int k=0; k<=8; ++k) {
            const double value=std::ldexp(static_cast<double>(p[k]),k);
            a[k]={value,value};
        }
        return interval_schur(a);
    };
    for (Int sign : {-1,1}) {
        Poly doubled{};
        doubled[8]=2;
        doubled[7]=4*sign;
        if (enclosed(doubled) != Decision::fallback ||
            exact_schur(doubled,g,false,0))
            throw std::runtime_error("strict equality missed");
        doubled[7]=2*sign;
        if (enclosed(doubled) != Decision::pass ||
            !exact_schur(doubled,g,false,0))
            throw std::runtime_error("strict interior missed");
        doubled[7]=6*sign;
        if (enclosed(doubled) != Decision::fail ||
            exact_schur(doubled,g,false,0))
            throw std::runtime_error("strict exterior missed");
    }
    volatile double tiny = std::numeric_limits<double>::denorm_min();
    const double twice = tiny + tiny;
    if (!(tiny > 0) || !(twice > tiny) || std::nextafter(0.,1.) != tiny)
        throw std::runtime_error("unexpected subnormal behavior");
    std::cout << "{\"status\":\"passed\",\"boundary_fallback\":2,\"strict_interior\":2,\"strict_exterior\":2,\"subnormal\":\"verified\"}"
              << '\n';
}
