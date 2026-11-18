"""Add deterministic host reductions of the verified terminal prefix path data."""
from pathlib import Path
import sys

source,target=map(Path,sys.argv[1:])
s=source.read_text()
def one(a,b):
    global s
    assert s.count(a)==1,a
    s=s.replace(a,b)

one('  let rbTransport: [[Double]]',
    '  let rbTransport: [[Double]]\n  let contact: [String: [UInt64]]')
one('rbTransport: [], dead: true','rbTransport: [], contact: [:], dead: true')
one('  let rbM=tau==2 ?',r'''
  // The top-two lift weights this rank/prefix statistic by L, rather than the
  // different first-merge reweighting used for functions of only the peak pair.
  var contact: [String: [UInt64]] = [:]
  if tau == 2 {
    precondition(n <= 32)
    let word = px.contents().bindMemory(to: UInt16.self,capacity:Np*1024)
    let rowMask = (UInt64(1)<<n)-1
    var open = [UInt64](repeating: rowMask,count:n)
    var queue = [Int](repeating: 0,count:N)
    var joint = [UInt64](repeating: 0,count:distanceM*distanceM)
    var ranks = [UInt64](repeating: 0,count:N)
    var lowGrows = UInt64(0)
    for particle in 0..<Np {
      for row in open.indices { open[row] = rowMask }
      let p=Int(word[particle*1024]), q=Int(word[particle*1024+1])
      let px=p%n, py=p/n, qx=q%n, qy=q/n
      precondition(p==Int(rb[4*particle+2]) && q==Int(rb[4*particle+3]))
      let sourceDist=abs(px-qx)+abs(py-qy)
      var j=0, firstNeighbor=0
      for k in 0..<(N-1) {
        if k==1 { continue }
        let v=Int(word[particle*1024+k]), x=v%n, y=v/n
        if abs(x-qx)+abs(y-qy)==1 { j=k; firstNeighbor=v; break }
        let bit=UInt64(1)<<x
        open[y] &= ~(bit | (bit>>1) | ((bit<<1)&rowMask))
        if y>0 { open[y-1] &= ~bit }
        if y+1<n { open[y+1] &= ~bit }
      }
      precondition(j>=2 && (open[qy] & (UInt64(1)<<qx)) != 0)
      let weight=UInt64(j-1)
      if (open[firstNeighbor/n] & (UInt64(1)<<(firstNeighbor%n))) != 0 {
        lowGrows += weight
      }
      open[qy] &= ~(UInt64(1)<<qx)
      queue[0]=q
      var head=0, tail=1, far=sourceDist
      while head<tail {
        let v=queue[head], x=v%n, y=v/n
        head+=1
        far=max(far,abs(px-x)+abs(py-y))
        let bit=UInt64(1)<<x
        if x>0 && (open[y] & (bit>>1)) != 0 {
          open[y] &= ~(bit>>1); queue[tail]=v-1; tail+=1
        }
        if x+1<n && (open[y] & (bit<<1)) != 0 {
          open[y] &= ~(bit<<1); queue[tail]=v+1; tail+=1
        }
        if y>0 && (open[y-1] & bit) != 0 {
          open[y-1] &= ~bit; queue[tail]=v-n; tail+=1
        }
        if y+1<n && (open[y+1] & bit) != 0 {
          open[y+1] &= ~bit; queue[tail]=v+n; tail+=1
        }
      }
      joint[sourceDist*distanceM+far]+=weight
      ranks[j]+=weight
    }
    precondition(ranks.reduce(0,+)==rawLiftL && joint.reduce(0,+)==rawLiftL)
    for d in 0..<distanceM {
      let own=joint[(d*distanceM)..<((d+1)*distanceM)].reduce(UInt64(0),+)
      let kernel=stride(from:d,to:groupDistanceHist.count,by:distanceM).reduce(UInt64(0)) { $0+groupDistanceHist[$1] }
      precondition(own==kernel)
    }
    contact=["contact_rank_L":ranks,"contact_joint_source_max_distance_L":joint,
             "first_contact_nonmerge_L":[lowGrows]]
  }
  let rbM=tau==2 ?''')
one('rbTransport: [nearest,fixed,adaptive,orientation].map{$0.map{$0/Double(Np)}}, dead: false',
    'rbTransport: [nearest,fixed,adaptive,orientation].map{$0.map{$0/Double(Np)}}, contact: contact, dead: false')
one('"transport_histograms":tau == 2 ? z.rbTransport : NSNull()]',
    '"transport_histograms":tau == 2 ? z.rbTransport : NSNull(), "contact_counts":z.contact]')
target.write_text(s)
