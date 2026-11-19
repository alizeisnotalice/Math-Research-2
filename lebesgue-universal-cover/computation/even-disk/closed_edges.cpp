#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>
using std::int32_t;
constexpr int cells = 480;
constexpr int velocities = 280;
constexpr int threshold = 460000000;
std::vector<int32_t> read_vector(const char* file) {
    std::vector<int32_t> values(cells * velocities);
    std::ifstream in(file, std::ios::binary);
    if (!in.read(reinterpret_cast<char*>(values.data()), values.size()*sizeof(int32_t)) || in.peek() != EOF) throw "incorrect vector size";
    return values;
}
int floor_half(int k) { return k >= 0 ? k/2 : -((-k+1)/2); }
int ceil_half(int k) { return -floor_half(-k); }
int main(int argc, char** argv) {
    if (argc != 3) return 1;
    const auto cost = read_vector(argv[1]), potential = read_vector(argv[2]);
    long long arcs = 0;
    for (int i=0; i<cells; ++i) for (int j=0; j<velocities; ++j) {
        int p=i*velocities+j;
        if (cost[p]<0 || potential[p]<0 || potential[p]>threshold) return 2;
        if (i<=145 && potential[p]>cost[p]/2) return 3;
        if (i>=334) { if (potential[p]!=threshold) return 4; continue; }
        for (int jj=0; jj<velocities; ++jj) {
            int first=std::max(0,ceil_half(2*i+j+jj-280)-1);
            int last=std::min(cells-1,floor_half(2*(i+1)+j+jj-278));
            for (int ii=first; ii<=last; ++ii) {
                int v=ii*velocities+jj;
                long long amount = ii>=334 ? cost[v]/2 : cost[v];
                if (potential[v]-static_cast<long long>(potential[p]) > amount) return 5;
                ++arcs;
            }
        }
    }
    std::cout << "states " << potential.size() << ", arcs " << arcs << ", threshold " << threshold << "/1000000000\n";
}
