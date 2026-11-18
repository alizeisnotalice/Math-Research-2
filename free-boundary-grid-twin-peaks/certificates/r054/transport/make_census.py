"""Derive short marked-path integer census from the r53 completion tables."""
from pathlib import Path
import sys

p,out=map(Path,sys.argv[1:])
s=p.read_text()
s=s.replace('// Exact census of the vertices reachable before the low peak\'s first neighbor.',
            '// Exact marked simple-path counts in early first-contact words.')
s=s.replace('using Key = std::tuple<int,int,int>;\nusing Row = std::array<Count,2>;\nusing Table = std::map<Key,Row>;',
            'using Key = std::tuple<int,int,int,int>;\nusing Table = std::map<Key,Count>;')
a=s.index('        std::vector<Row> weights(size);')
b=s.index('        return out;',a)
s=s[:a]+'''        std::vector<Count> weights(size), ways(size);
        Table out;
        const std::set<int> thresholds{std::max(3,n_/2+1),n_-1};
        for (int p=0;p<static_cast<int>(verts.size());++p) {
            const int P=verts[p], d=distance(P,q);
            if (d>=*thresholds.rbegin()) continue;
            std::fill(ways.begin(),ways.end(),0);
            ways[bit(p)]=1;
            for (Mask s=bit(p);s<size;++s) {
                if (!ways[s]) continue;
                const int len=std::popcount(s);
                if (len<6) {
                    for (Mask v=fringe[s]&~s;v;drop(v)) ways[s|bit(first(v))]+=ways[s];
                }
                const Mask K=expanded[s];
                const Mask U=all_&~(K|ring_[K]);
                if (!weights[s]) {
                    for (Mask v=edges_[q];v;drop(v)) weights[s]+=suffix_[K|bit(q)|bit(first(v))];
                    assert(weights[s]>0);
                }
                const Count coefficient=ways[s]*len*orbit*weights[s];
                for (int t:thresholds) {
                    if (d>=t) continue;
                    std::array<std::uint64_t,5> counts{};
                    Mask choices=U;
                    for (Mask a=U;a;drop(a)) {
                        const int v=first(a);
                        if (distance(P,v)<t && (distance(q,v)%2)==0) choices&=~bit(v);
                    }
                    paths(q,choices,P,t,0,counts);
                    for (int k=1;k<=4;++k) {
                        if (counts[k]) out[{len+1,d,t,k}]+=coefficient*counts[k];
                    }
                }
            }
        }
''' +s[b:]
a=s.index('    int distance(')
s=s[:a]+'''    void paths(int at,Mask remain,int p,int threshold,int length,
               std::array<std::uint64_t,5>& counts) const {
        if (length && distance(p,at)>=threshold) ++counts[length];
        if (length==4) return;
        for (Mask v=edges_[at]&remain;v;drop(v)) {
            const int u=first(v);
            paths(u,remain&~bit(u),p,threshold,length+1,counts);
        }
    }

'''+s[a:]
s=s.replace('                    for (int a=0;a<2;++a) partial[i][key][a]+=row[a];','                    partial[i][key]+=row;')
s=s.replace('            for (int a=0;a<2;++a) rows[key][a]+=row[a];','            rows[key]+=row;')
s=s.replace('const auto [rank,d,far]=key;', 'const auto [rank,d,t,k]=key;')
s=s.replace("<<far<<' '<<decimal(row[0])<<' '<<decimal(row[1])", "<<t<<' '<<k<<' '<<decimal(row)")
out.write_text(s)
