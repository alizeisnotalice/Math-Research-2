#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"
#define main search_program_main
#include "psi_search_undo.cpp"
#undef main
#pragma clang diagnostic pop
#include <numeric>
static bool arrays_equal(const State& a,const State& b){return a.a==b.a&&a.pos==b.pos&&a.l==b.l&&a.b==b.b&&a.sole==b.sole;}
int main(){
 for(int n:{2,3}){
  int N=n*n;std::vector<int>w(N);std::iota(w.begin(),w.end(),0);
  uint64_t words=0,valid=0,invalid=0,noops=0;
  do{
   State s(n,w[0]);s.a=w;for(int j=0;j<N;++j)s.pos[w[j]]=j+1;
   bool one=true;for(int v=0;v<N;++v)if(v!=s.root&&s.earlycount(v)==0)one=false;
   if(!one)continue;s.rebuild();++words;
   for(int i=0;i<N;++i)for(int j=0;j<N;++j){
    State x=s;State::Undo u;bool ok=x.swap_try(i,j,u);
    if(ok){++valid;State rebuilt=x;rebuilt.rebuild();if(!arrays_equal(x,rebuilt)||std::abs(x.psi-rebuilt.psi)>1e-11*std::max(1.,std::abs(rebuilt.psi)))return 2;x=std::move(rebuilt);}
    else if(u.i<0)++noops;else ++invalid;
    x.rollback(u);
    if(!arrays_equal(x,s)||x.psi!=s.psi)return 3;
   }
  }while(std::next_permutation(w.begin(),w.end()));
  std::cout<<"n="<<n<<" onepeak_words="<<words<<" valid_candidates="<<valid<<" invalid_candidates="<<invalid<<" noops="<<noops<<" rollback_checks="<<(valid+invalid+noops)<<" PASS\n";
 }
}
