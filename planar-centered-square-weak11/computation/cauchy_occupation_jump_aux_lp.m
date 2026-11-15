function cauchy_occupation_jump_aux_lp(h,delta,N,subdiv,output_root)
% R62 V3: telescope adjacent density jumps before interval curvature bounds.
% PROPOSAL ONLY: nondirected MATLAB. Root owns batch lifecycle.
assert(h>0&&delta>0&&N>=1&&subdiv>=1&&subdiv==floor(subdiv));
assert(startsWith(output_root,'/Volumes/data/'));if ~isfolder(output_root),mkdir(output_root);end
started=tic;orbits=struct('center',{},'points',{});rep=[];J=0;m=2*N+1;orbmap=zeros(m,m);
for i=0:N,for j=0:i
 z=delta*[i,j];pts=[];for sx=[-1,1],for sy=[-1,1],pts=[pts;sx*z(1),sy*z(2);sx*z(2),sy*z(1)];end,end %#ok<AGROW>
 J=J+1;orbits(J).center=z;orbits(J).points=unique(pts,'rows');rep=[rep;z]; %#ok<AGROW>
 for p=1:size(orbits(J).points,1),u=round(orbits(J).points(p,:)/delta)+N+1;orbmap(u(1),u(2))=J;end
end,end
assert(all(orbmap(:)>0));sizes=arrayfun(@(o)size(o.points,1),orbits)';
% Vertical and horizontal physical interfaces. Exterior densities are zero.
Ex=m*(m+1);Ey=Ex;Etot=Ex+Ey;xe=zeros(Ex,1);xy=xe;ye=xe;yx=xe;
JD=spalloc(Etot,J+1,2*Etot);r=0;
for k=1:m+1,for t=1:m
 r=r+1;xe(r)=delta*(k-N-1.5);xy(r)=delta*(t-N-1);
 if k<=m,JD(r,1+orbmap(k,t))=JD(r,1+orbmap(k,t))+1;end
 if k>=2,JD(r,1+orbmap(k-1,t))=JD(r,1+orbmap(k-1,t))-1;end
 ye(r)=xe(r);yx(r)=xy(r);
 if k<=m,JD(Ex+r,1+orbmap(t,k))=JD(Ex+r,1+orbmap(t,k))+1;end
 if k>=2,JD(Ex+r,1+orbmap(t,k-1))=JD(Ex+r,1+orbmap(t,k-1))-1;end
end,end
Bcount=J*subdiv^2;R=4*Bcount;V=zeros(R,J+1);target=zeros(R,1);ox=zeros(Bcount,1);oy=ox;
PxL=zeros(Bcount,Ex);PxU=PxL;PyL=zeros(Bcount,Ey);PyU=PyL;bi=0;
for i=1:J
 for ix=0:subdiv-1,for iy=0:subdiv-1
  bi=bi+1;ids=4*(bi-1)+(1:4);
  B=[rep(i,1)-delta/2+delta*[ix,ix+1]/subdiv;rep(i,2)-delta/2+delta*[iy,iy+1]/subdiv];
  vertices=[B(1,1),B(2,1);B(1,1),B(2,2);B(1,2),B(2,1);B(1,2),B(2,2)];
  V(ids,1)=ph(vertices(:,1),vertices(:,2),h);target(ids)=i;
  [px,pxx]=cauchyjet(B(1,:),h);[py,pyy]=cauchyjet(B(2,:),h);
  ux=imul(pxx,py);uy=imul(px,pyy);ox(bi)=ux(2);oy(bi)=uy(2);
  for j=1:J
   pts=orbits(j).points;
   for k=1:size(pts,1)
    V(ids,j+1)=V(ids,j+1)+boxprob(abs(vertices(:,1)-pts(k,1)),delta,h).*boxprob(abs(vertices(:,2)-pts(k,2)),delta,h);
   end
  end
  [kl,ku]=kprimevec(B(1,1)-xe,B(1,2)-xe,h);
  [tl,tu]=boxrangevec(B(2,1)-xy,B(2,2)-xy,delta,h);
  [pl,pu]=mulvec(kl,ku,tl,tu);PxL(bi,:)=pl';PxU(bi,:)=pu';
  [kl,ku]=kprimevec(B(2,1)-ye,B(2,2)-ye,h);
  [tl,tu]=boxrangevec(B(1,1)-yx,B(1,2)-yx,delta,h);
  [pl,pu]=mulvec(kl,ku,tl,tu);PyL(bi,:)=pl';PyU(bi,:)=pu';
 end,end
end
source_margin=1e-7;V(:,1)=V(:,1)-source_margin*max(1,abs(V(:,1)));
baseA=-V;baseA(sub2ind(size(baseA),(1:R)',target+1))=baseA(sub2ind(size(baseA),(1:R)',target+1))+1;
ell=delta/subdiv;E=sparse((1:R)',repelem((1:Bcount)',4),repmat(ell^2/8,R,1),R,Bcount);
% Variable order: base(lambda,a), dplus(allinterfaces), dminus, eta_x,eta_y.
sourceA=[sparse(baseA),sparse(R,2*Etot),E,E];
Cx=sparse((1:Bcount)',ones(Bcount,1),ox,Bcount,J+1);Cy=sparse((1:Bcount)',ones(Bcount,1),oy,Bcount,J+1);
ZX=sparse(Bcount,Ex);ZY=sparse(Bcount,Ey);ZZ=sparse(Bcount,Bcount);II=speye(Bcount);
curvA=[Cx,sparse(PxU),ZY,sparse(-PxL),ZY,-II,ZZ;Cy,ZX,sparse(PyU),ZX,sparse(-PyL),ZZ,-II];
Aeq=[JD,-speye(Etot),speye(Etot),sparse(Etot,2*Bcount)];beq=zeros(Etot,1);
[gn,gw]=gauss(24);[x,y]=diamond_points(24,96);D=kmatrix(x,y,h,delta,orbits,gn,gw);q=shortQ(x,y,h);
Amat=[sourceA;curvA;sparse([-q,-D]),sparse(numel(x),2*Etot+2*Bcount)];bv=[zeros(R+2*Bcount,1);-ones(numel(x),1)];
cost=pi*h*[1;delta^2*sizes;zeros(2*Etot+2*Bcount,1)];
opts=optimoptions('linprog','Display','none','ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'MaxIterations',20000);
[v,lpval,ef,lpout]=linprog(cost,Amat,bv,Aeq,beq,zeros(numel(cost),1),[],opts);assert(ef>0,'jump auxiliary LP failed');v=max(v,0);
[fx,fy]=diamond_points(48,192);FD=kmatrix(fx,fy,h,delta,orbits,gn,gw);fv=v(1)*shortQ(fx,fy,h)+FD*v(2:J+1);[dm,idx]=min(fv);assert(dm>0);v=v/dm;
actual_jump=JD*v(1:J+1);jump_parts=v(J+2:J+1+Etot)+v(J+2+Etot:J+1+2*Etot);
out=struct('status','completed_numeric_telescoped_jump_aux_PROPOSAL','matlab_version',version,'threads',1,'h',h,'delta',delta,'N',N,'subdiv',subdiv,'source_margin',source_margin,'domain_halfwidth',(N+.5)*delta,'D4_orbit_count',J,'physical_cells',sum(sizes),'source_subrectangles',Bcount,'source_rows',R,'curvature_rows',2*Bcount,'jump_equality_rows',Etot,'LP_variables',numel(v),'LP_nnz',nnz(Amat)+nnz(Aeq),'source_origin_mass',v(1),'density_coefficients',v(2:J+1),'orbit_representatives',rep,'orbit_sizes',sizes,'LP_objective_before_D_rescale',lpval,'LP_exitflag',ef,'LP_output',lpout,'D_coarse_count',numel(x),'D_fine_count',numel(fx),'D_min_before_rescale',dm,'D_minimum_point',[fx(idx),fy(idx)],'cost_from_exact_mass_formula_float',cost'*v,'source_constraint_min_float',-max(sourceA*v),'curvature_constraint_min_float',-max(curvA*v),'jump_equality_max_error_float',max(abs(Aeq*v)),'density_jump_l1_float',sum(abs(actual_jump)),'jump_part_extra_float',sum(jump_parts-abs(actual_jump)),'Qcell_gauss_order',24,'elapsed_seconds',toc(started),'scope','Telescope adjacent physical-cell density jumps before curvature intervals. Floating proposal only; no general upper without exact source/full-D certificates.');
fid=fopen(fullfile(output_root,'occupation-jump-aux-lp.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
fprintf('JUMP-AUX h=%.6g delta=%.6g N=%d s=%d J=%d interfaces=%d vars=%d cost=%.15g source=%.8g curvature=%.8g equality=%.5g elapsed=%.4g\n',h,delta,N,subdiv,J,Etot,numel(v),cost'*v,out.source_constraint_min_float,out.curvature_constraint_min_float,out.jump_equality_max_error_float,toc(started));
end
function [p,dd]=cauchyjet(x,h)
dmin=min(abs(x));if x(1)<=0&&x(2)>=0,dmin=0;end;dmax=max(abs(x));
p=[h/(pi*(h*h+dmax*dmax)),h/(pi*(h*h+dmin*dmin))];
z=[x(1),x(2)];for c=[-h,0,h],if x(1)<=c&&c<=x(2),z(end+1)=c;end,end %#ok<AGROW>
v=2*h*(3*z.^2-h*h)./(pi*(h*h+z.^2).^3);dd=[min(v),max(v)];
end
function [lo,hi]=kprimevec(xl,xr,h)
a=-2*h*xl./(pi*(h*h+xl.^2).^2);b=-2*h*xr./(pi*(h*h+xr.^2).^2);lo=min(a,b);hi=max(a,b);
c=h/sqrt(3);mag=9/(8*sqrt(3)*pi*h*h);mask=xl<=c&c<=xr;lo(mask)=-mag;mask=xl<=-c&-c<=xr;hi(mask)=mag;
end
function [lo,hi]=boxrangevec(xl,xr,d,h)
dmin=min(abs(xl),abs(xr));dmin(xl<=0&xr>=0)=0;dmax=max(abs(xl),abs(xr));lo=boxprob(dmax,d,h);hi=boxprob(dmin,d,h);
end
function [lo,hi]=mulvec(al,au,bl,bu)
z=[al.*bl,al.*bu,au.*bl,au.*bu];lo=min(z,[],2);hi=max(z,[],2);
end
function y=imul(a,b)
z=[a(1)*b(1),a(1)*b(2),a(2)*b(1),a(2)*b(2)];y=[min(z),max(z)];
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
