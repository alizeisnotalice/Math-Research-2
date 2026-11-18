#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
using I=long long;
struct Box{std::array<I,11> a,b;};
struct Counts{I total=0,keep=0,combo=0;};
std::vector<unsigned char> contents; size_t pos=0;
constexpr I unit=I(1)<<54;
const int width[4]={700,670,606,575};const int ax[4]={0,3,6,9};

bool outside(const Box& t){
 if((__int128)t.a[2]*7*15>(__int128)22*unit)return true;
 for(int i=0;i<4;i++){
   int d=ax[i];I x=std::max({t.a[d],-t.b[d],I(0)}),y=std::max({t.a[d+1],-t.b[d+1],I(0)});
   // widths use 1/10000. Compare only after outward-safe dyadic rounding to 2^-30.
   x>>=24;y>>=24;I scale=I(1)<<30;
   if((__int128(1)*x*x+__int128(1)*y*y)*100000000>(__int128)width[i]*width[i]*scale*scale)return true;
 }
 return false;
}
void walk(Box& t,Counts &s,bool out){
 unsigned code=contents.at(pos++);if(code>=253){s.total++;if(!out&&!outside(t)){s.keep++;if(code==253)s.combo++;}return;}
 assert(code<11);I a=t.a[code],b=t.b[code];assert(((a+b)&1)==0);I mid=(a+b)/2;
 t.b[code]=mid;walk(t,s,out||outside(t));t.b[code]=b;t.a[code]=mid;walk(t,s,out||outside(t));t.a[code]=a;
}
int main(int argc,char** argv){assert(argc==24);Box root;for(int i=0;i<11;i++){root.a[i]=std::stoll(argv[2+i]);root.b[i]=std::stoll(argv[13+i]);}Counts sum;
 for(int part=0;part<2048;part++){char file[32];sprintf(file,"/part%04d.npy",part);std::ifstream f(std::string(argv[1])+file,std::ios::binary);assert(f);contents.assign(std::istreambuf_iterator<char>(f),{});assert(contents.size()>=10&&contents[0]==147&&contents[6]==1);pos=10+contents[8]+256*contents[9];std::string header(contents.begin()+10,contents.begin()+pos);assert(header.find("|u1")!=std::string::npos&&header.find("False")!=std::string::npos);Box t=root;
  for(int i=0;i<11;i++){I m=(t.a[i]+t.b[i])/2;if((part>>(10-i))&1)t.a[i]=m;else t.b[i]=m;}
  walk(t,sum,outside(t));assert(pos==contents.size());}
 printf("{\"leaves\":%lld,\"retained\":%lld,\"retained_multi\":%lld}\n",sum.total,sum.keep,sum.combo);
}
