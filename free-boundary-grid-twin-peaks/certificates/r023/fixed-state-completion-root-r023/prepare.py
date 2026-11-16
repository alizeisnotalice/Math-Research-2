from pathlib import Path
import hashlib,json
src=Path('/Volumes/data/CodexTemp/T2/a022-main-rounds/free-boundary-grid-twin-peaks/certificates/r021-root-n7/root_counts_d4.cpp')
s=src.read_text();oldhash=hashlib.sha256(src.read_bytes()).hexdigest()
s=s.replace('static int classes=0;','static int classes=0;\nconstexpr int MAXSEEDS=8;\nstatic std::vector<Mask> seeds, bans;\nstatic std::vector<std::string> seed_names;')
s=s.replace('  U256 sq{},fp{},fm{},gp{},gm{};','  U256 sq{},fp{},fm{},gp{},gm{};\n  std::array<U256,MAXSEEDS> seedp{},seedm{};')
s=s.replace('  for(int v=0;v<N;v++)assert(vclass[v]<classes);','''  for(int v=0;v<N;v++)assert(vclass[v]<classes);
  auto seed=[&](Mask m,std::string name){seeds.push_back(m);bans.push_back(neighbourhood(m));seed_names.push_back(name);};
  seed(Mask(1),"corner-root");seed(Mask(1)<<(n+1),"inner-root");
  for(int a=1;a<=std::min(3,n-2);a++){
    Mask m=1;for(int y=1;y<=a;y++)for(int x=1;x<=a;x++)m|=Mask(1)<<(y*n+x);
    seed(m,"block"+std::to_string(a)+"-dormant");seed(m|2,"block"+std::to_string(a)+"-bridge");
  }
  assert(seeds.size()<=MAXSEEDS);''')
s=s.replace('  U256 Sp{},Sm{},Fp{},Fm{},Gp{},Gm{};','  U256 Sp{},Sm{},Fp{},Fm{},Gp{},Gm{};\n  std::array<U256,MAXSEEDS> SeedP{},SeedM{};for(size_t j=0;j<seeds.size();j++)SeedP[j]=fact;')
needle='          U256 w=mulsmall(cur[z],layers[k].orb[z]); add(acc[id].t[k],w);'
insert='''
          // Every distinct orbit member occurs 8/orbit_size times among
          // images[0..7]. Multiply BEFORE division to retain stabilizers.
          for(size_t j=0;j<seeds.size();j++){
            unsigned hits=0;for(int op=0;op<8;op++)hits+=!(images[op]&bans[j]);
            unsigned mult=hits*layers[k].orb[z];assert(mult%8==0);mult/=8;
            if(mult){U256 term=mulsmall(cur[z],mult);if(k&1)add(acc[id].seedm[j],term);else add(acc[id].seedp[j],term);}
          }
'''
assert s.count(needle)==1;s=s.replace(needle,needle+insert)
needle='      add(TkP[k],acc[i].t[k]);';assert s.count(needle)==1;s=s.replace(needle,needle+'\n      for(size_t j=0;j<seeds.size();j++){add(SeedP[j],acc[i].seedp[j]);add(SeedM[j],acc[i].seedm[j]);}')
needle='  std::cout<<"{\\\"n\\\":"<<n'
assert s.count(needle)==1
insert='''  std::cout<<"{\\\"fixed_state_completions\\\":[";
  for(size_t j=0;j<seeds.size();j++){
    assert(cmp(SeedP[j],SeedM[j])>=0);U256 h=sub(SeedP[j],SeedM[j]);
    unsigned sz=__builtin_popcountll(seeds[j]);for(unsigned t=0;t<sz;t++)h=divsmall(h,N-t);
    if(j==0)assert(str(h)==sstr(RootP[0],RootM[0]));
    if(j==1)assert(str(h)==sstr(RootP[n+1],RootM[n+1]));
    if(j)std::cout<<',';
    std::cout<<"{\\\"name\\\":\\\""<<seed_names[j]<<"\\\",\\\"S\\\":"<<seeds[j]<<",\\\"size\\\":"<<sz<<",\\\"ban\\\":"<<bans[j]<<",\\\"completion_count\\\":\\\""<<str(h)<<"\\\"}";
  }
  std::cout<<"],";
'''
s=s.replace(needle,insert+needle.replace('"{\\\"n','"\\\"n'))
Path('fixed_state_d4.cpp').write_text(s)
Path('source-provenance.json').write_text(json.dumps({'base_source':str(src),'base_sha256':oldhash,'new_sha256':hashlib.sha256(s.encode()).hexdigest(),'change':'new seed forbidden-orbit counters only; existing D4 values/recurrence unchanged'},indent=2)+'\n')
