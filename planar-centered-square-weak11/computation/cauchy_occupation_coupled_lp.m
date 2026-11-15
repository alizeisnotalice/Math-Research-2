function cauchy_occupation_coupled_lp(h,delta,N,subdiv,robust,output_root)
% R62 coupled whole-cell source LP PROPOSAL. Nondirected MATLAB arithmetic.
% Root launches/owns batch processes. Exact postchecker is separate.
assert(h>0&&delta>0&&N>=1&&subdiv>=1&&subdiv==floor(subdiv));
assert(ismember(robust,[0,1]));assert(startsWith(output_root,'/Volumes/data/'));
if ~isfolder(output_root),mkdir(output_root);end
started=tic;orbits=struct('center',{},'points',{});rep=[];J=0;
for i=0:N,for j=0:i
 z=delta*[i,j];pts=[];for sx=[-1,1],for sy=[-1,1],pts=[pts;sx*z(1),sy*z(2);sx*z(2),sy*z(1)];end,end %#ok<AGROW>
 J=J+1;orbits(J).center=z;orbits(J).points=unique(pts,'rows');rep=[rep;z]; %#ok<AGROW>
end,end
sizes=arrayfun(@(o)size(o.points,1),orbits)';
% Rows retain full linear combination at shared corners. Within each B,
% bilinear interpolation error <= (ell_x^2 Mxx+ell_y^2 Myy)/8.
R=4*J*subdiv^2;S=zeros(R,J+1);target=zeros(R,1);row=0;
for i=1:J
 for ix=0:subdiv-1,for iy=0:subdiv-1
  B=[rep(i,1)-delta/2+delta*[ix,ix+1]/subdiv;rep(i,2)-delta/2+delta*[iy,iy+1]/subdiv];
  vertices=[B(1,1),B(2,1);B(1,1),B(2,2);B(1,2),B(2,1);B(1,2),B(2,2)];
  val=zeros(4,J+1);val(:,1)=ph(vertices(:,1),vertices(:,2),h);
  E=zeros(1,J+1);
  if robust
   [px,pxx]=cauchyjet(B(1,:),h);[py,pyy]=cauchyjet(B(2,:),h);
   ux=imul(pxx,py);uy=imul(px,pyy);
   E(1)=delta^2/(8*subdiv^2)*(max(0,ux(2))+max(0,uy(2)));
  end
  for j=1:J
   pts=orbits(j).points;ux=[0,0];uy=[0,0];
   for k=1:size(pts,1)
    val(:,j+1)=val(:,j+1)+boxprob(abs(vertices(:,1)-pts(k,1)),delta,h).*boxprob(abs(vertices(:,2)-pts(k,2)),delta,h);
    if robust
     [px,pxx]=boxjet(B(1,:)-pts(k,1),delta,h);[py,pyy]=boxjet(B(2,:)-pts(k,2),delta,h);
     ux=ux+imul(pxx,py);uy=uy+imul(px,pyy);
    end
   end
   if robust,E(j+1)=delta^2/(8*subdiv^2)*(max(0,ux(2))+max(0,uy(2)));end
  end
  ids=row+(1:4);S(ids,:)=val-E;target(ids)=i;row=row+4;
 end,end
end
% Mild origin-variable slack makes exact postverification/repair practical.
source_margin=1e-7;S(:,1)=S(:,1)-source_margin*max(1,abs(S(:,1)));
sourceA=-S;sourceA(sub2ind(size(sourceA),(1:R)',target+1))=sourceA(sub2ind(size(sourceA),(1:R)',target+1))+1;
[gn,gw]=gauss(24);[x,y]=diamond_points(24,96);D=kmatrix(x,y,h,delta,orbits,gn,gw);q=shortQ(x,y,h);
Amat=[sourceA;-q,-D];bv=[zeros(R,1);-ones(numel(x),1)];cost=pi*h*[1;delta^2*sizes];
opts=optimoptions('linprog','Display','none','ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'MaxIterations',20000);
[v,lpval,ef,lpout]=linprog(cost,Amat,bv,[],[],zeros(J+1,1),[],opts);
assert(ef>0,'coupled occupation LP failed');v=max(v,0);
[fx,fy]=diamond_points(48,192);FD=kmatrix(fx,fy,h,delta,orbits,gn,gw);fv=v(1)*shortQ(fx,fy,h)+FD*v(2:end);[dm,idx]=min(fv);assert(dm>0);v=v/dm;
out=struct('status','completed_numeric_coupled_source_PROPOSAL','matlab_version',version,'threads',1,'h',h,'delta',delta,'N',N,'subdiv',subdiv,'robust_curvature_rows',robust,'source_margin',source_margin,'domain_halfwidth',(N+.5)*delta,'D4_orbit_count',J,'physical_cells',sum(sizes),'source_rows',R,'source_origin_mass',v(1),'density_coefficients',v(2:end),'orbit_representatives',rep,'orbit_sizes',sizes,'LP_objective_before_D_rescale',lpval,'LP_exitflag',ef,'LP_output',lpout,'D_coarse_count',numel(x),'D_fine_count',numel(fx),'D_min_before_rescale',dm,'D_minimum_point',[fx(idx),fy(idx)],'cost_from_exact_mass_formula_float',cost'*v,'source_constraint_min_float',min(S*v-v(target+1)),'Qcell_gauss_order',24,'elapsed_seconds',toc(started),'scope','Whole-cell bilinear-interpolation source constraints; floating proposal only. No new general upper without exact source and full-D certificates.');
fid=fopen(fullfile(output_root,'occupation-coupled-lp.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
fprintf('COUPLED h=%.6g delta=%.6g N=%d s=%d robust=%d J=%d rows=%d cost=%.15g residual=%.8g elapsed=%.4g\n',h,delta,N,subdiv,robust,J,R,cost'*v,out.source_constraint_min_float,toc(started));
end
function [p,dd]=cauchyjet(x,h)
sq=squareI(x);den=h*h+sq;p=h/pi./fliplr(den);
num=2*h/pi*(3*sq-h*h);dd=idiv(num,den.^3);
end
function [p,dd]=boxjet(x,d,h)
sq=squareI(x);p=[boxprob(sqrt(sq(2)),d,h),boxprob(sqrt(sq(1)),d,h)];
a=kprime(x+d/2,h);b=kprime(x-d/2,h);dd=[a(1)-b(2),a(2)-b(1)];
end
function a=kprime(x,h)
den=(h*h+squareI(x)).^2;a=idiv(-2*h/pi*fliplr(x),den);
end
function y=squareI(x)
y=[min(x.^2),max(x.^2)];if x(1)<=0&&x(2)>=0,y(1)=0;end
end
function y=imul(a,b)
z=[a(1)*b(1),a(1)*b(2),a(2)*b(1),a(2)*b(2)];y=[min(z),max(z)];
end
function y=idiv(a,b)
assert(b(1)>0);y=imul(a,1./fliplr(b));
end
function p=ph(x,y,h)
p=h*h./(pi^2*(h*h+x.^2).*(h*h+y.^2));
end
function p=boxprob(dist,delta,s)
p=atan2(s*delta,s*s+dist.^2-delta^2/4)/pi;
end
function D=kmatrix(x,y,h,delta,orbits,gn,gw)
D=zeros(numel(x),numel(orbits));
for j=1:numel(orbits)
 pts=orbits(j).points;v=zeros(size(x));
 for k=1:size(pts,1)
  dx=abs(x-pts(k,1));dy=abs(y-pts(k,2));
  for q=1:numel(gn),s=h*gn(q);v=v+2*pi*h*gw(q)*boxprob(dx,delta,s).*boxprob(dy,delta,s);end
 end
 D(:,j)=v;
end
end
function q=shortQ(a,b,t)
a=abs(a);b=abs(b);lo=min(a,b);hi=max(a,b);h=hi-lo;den=t*t+lo.*hi;z=t*h./den;ar=ones(size(z));ii=z~=0;ar(ii)=atan(z(ii))./z(ii);q=(2/pi)*(atan(t./hi)-lo*t./den.*ar)./(lo+hi);
ii=lo>8*t;
if any(ii(:)),[nodes,wt]=gauss(16);aa=lo(ii).^2;bb=hi(ii).^2;v=zeros(size(aa));for j=1:numel(nodes),s2=(t*nodes(j))^2;v=v+wt(j)*s2./((s2+aa).*(s2+bb));end;q(ii)=2*t/pi*v;end
end
function [x,w]=gauss(n)
b=(1:n-1)./sqrt(4*(1:n-1).^2-1);[V,D]=eig(diag(b,1)+diag(b,-1));[x,ix]=sort(diag(D));x=(x+1)/2;w=V(1,ix).^2;
end
function [x,y]=diamond_points(N,Nb)
x=[];y=[];for i=0:N,for j=0:i,if i+j<=N&&i+j>0,x(end+1,1)=i/N;y(end+1,1)=j/N;end,end,end %#ok<AGROW>
h=(0:Nb)'/(2*Nb);p=unique([x,y;1-h,h],'rows');x=p(:,1);y=p(:,2);
end
