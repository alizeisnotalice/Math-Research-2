function cauchy_occupation_tensor_jump_lp(h,grid_input,subdiv,output_root)
% R65: nonuniform symmetric tensor rectangles, exact boundary input strings.
% Ordinary double LP PROPOSAL only. Root owns launch/resource/exit control.
assert(h>0&&subdiv>=1&&subdiv==floor(subdiv));assert(startsWith(output_root,'/Volumes/data/'));
if ~isfolder(output_root),mkdir(output_root);end
started=tic;raw=jsondecode(fileread(grid_input));tokens=string(raw.axis_edges);ed=zeros(numel(tokens),1);
for k=1:numel(tokens),ed(k)=fraction_value(tokens(k));end
m=numel(ed)-1;assert(mod(m,2)==1&&all(diff(ed)>0)&&max(abs(ed+flipud(ed)))<1e-14);
width=diff(ed);center=(ed(1:end-1)+ed(2:end))/2;mid=(m+1)/2;assert(center(mid)==0);
J=0;orbits=struct('center',{},'width',{},'cells',{});rep=[];repwidth=[];orbmap=zeros(m,m);
for i=mid:m,for j=mid:i
 ij=[];
 for sx=[-1,1],for sy=[-1,1]
  if sx==1,ii=i;else,ii=m+1-i;end;if sy==1,jj=j;else,jj=m+1-j;end
  ij=[ij;ii,jj;jj,ii]; %#ok<AGROW>
 end,end
 ij=unique(ij,'rows');cells=[center(ij(:,1)),center(ij(:,2)),width(ij(:,1)),width(ij(:,2))];
 J=J+1;orbits(J).center=[center(i),center(j)];orbits(J).width=[width(i),width(j)];orbits(J).cells=cells;
 rep=[rep;orbits(J).center];repwidth=[repwidth;orbits(J).width]; %#ok<AGROW>
 for p=1:size(ij,1),assert(orbmap(ij(p,1),ij(p,2))==0);orbmap(ij(p,1),ij(p,2))=J;end
end,end
assert(all(orbmap(:)>0));sizes=arrayfun(@(o)size(o.cells,1),orbits)';areas=repwidth(:,1).*repwidth(:,2);
% Each vertical interface shares a source row with its ACTUAL y width.
Ex=m*(m+1);Ey=Ex;Etot=Ex+Ey;xe=zeros(Ex,1);xy=xe;xyw=xe;ye=xe;yx=xe;yxw=xe;
JD=spalloc(Etot,J+1,2*Etot);r=0;
for k=1:m+1,for t=1:m
 r=r+1;xe(r)=ed(k);xy(r)=center(t);xyw(r)=width(t);
 if k<=m,JD(r,1+orbmap(k,t))=JD(r,1+orbmap(k,t))+1;end
 if k>=2,JD(r,1+orbmap(k-1,t))=JD(r,1+orbmap(k-1,t))-1;end
 ye(r)=ed(k);yx(r)=center(t);yxw(r)=width(t);
 if k<=m,JD(Ex+r,1+orbmap(t,k))=JD(Ex+r,1+orbmap(t,k))+1;end
 if k>=2,JD(Ex+r,1+orbmap(t,k-1))=JD(Ex+r,1+orbmap(t,k-1))-1;end
end,end
Bcount=J*subdiv^2;R=4*Bcount;V=zeros(R,J+1);target=zeros(R,1);ox=zeros(Bcount,1);oy=ox;
PxL=zeros(Bcount,Ex);PxU=PxL;PyL=zeros(Bcount,Ey);PyU=PyL;errx=zeros(R,1);erry=errx;bi=0;
for i=1:J
 wx=repwidth(i,1);wy=repwidth(i,2);
 for ix=0:subdiv-1,for iy=0:subdiv-1
  bi=bi+1;ids=4*(bi-1)+(1:4);
  B=[rep(i,1)-wx/2+wx*[ix,ix+1]/subdiv;rep(i,2)-wy/2+wy*[iy,iy+1]/subdiv];
  vertices=[B(1,1),B(2,1);B(1,1),B(2,2);B(1,2),B(2,1);B(1,2),B(2,2)];
  V(ids,1)=ph(vertices(:,1),vertices(:,2),h);target(ids)=i;
  [px,pxx]=cauchyjet(B(1,:),h);[py,pyy]=cauchyjet(B(2,:),h);ux=imul(pxx,py);uy=imul(px,pyy);ox(bi)=ux(2);oy(bi)=uy(2);
  errx(ids)=(wx/subdiv)^2/8;erry(ids)=(wy/subdiv)^2/8;
  for j=1:J
   cells=orbits(j).cells;
   for k=1:size(cells,1)
    V(ids,j+1)=V(ids,j+1)+boxprob(abs(vertices(:,1)-cells(k,1)),cells(k,3),h).*boxprob(abs(vertices(:,2)-cells(k,2)),cells(k,4),h);
   end
  end
  [kl,ku]=kprimevec(B(1,1)-xe,B(1,2)-xe,h);[tl,tu]=boxrangevec(B(2,1)-xy,B(2,2)-xy,xyw,h);
  [pl,pu]=mulvec(kl,ku,tl,tu);PxL(bi,:)=pl';PxU(bi,:)=pu';
  [kl,ku]=kprimevec(B(2,1)-ye,B(2,2)-ye,h);[tl,tu]=boxrangevec(B(1,1)-yx,B(1,2)-yx,yxw,h);
  [pl,pu]=mulvec(kl,ku,tl,tu);PyL(bi,:)=pl';PyU(bi,:)=pu';
 end,end
end
source_margin=1e-7;V(:,1)=V(:,1)-source_margin*max(1,abs(V(:,1)));
baseA=-V;baseA(sub2ind(size(baseA),(1:R)',target+1))=baseA(sub2ind(size(baseA),(1:R)',target+1))+1;
% Anisotropic interpolation errors use each target subrectangle's own lengths.
EX=sparse((1:R)',repelem((1:Bcount)',4),errx,R,Bcount);EY=sparse((1:R)',repelem((1:Bcount)',4),erry,R,Bcount);
sourceA=[sparse(baseA),sparse(R,2*Etot),EX,EY];
Cx=sparse((1:Bcount)',ones(Bcount,1),ox,Bcount,J+1);Cy=sparse((1:Bcount)',ones(Bcount,1),oy,Bcount,J+1);
ZX=sparse(Bcount,Ex);ZY=sparse(Bcount,Ey);ZZ=sparse(Bcount,Bcount);II=speye(Bcount);
curvA=[Cx,sparse(PxU),ZY,sparse(-PxL),ZY,-II,ZZ;Cy,ZX,sparse(PyU),ZX,sparse(-PyL),ZZ,-II];
Aeq=[JD,-speye(Etot),speye(Etot),sparse(Etot,2*Bcount)];beq=zeros(Etot,1);
[gn,gw]=gauss(24);[x,y]=diamond_points(24,96);D=kmatrix(x,y,h,orbits,gn,gw);q=shortQ(x,y,h);
Amat=[sourceA;curvA;sparse([-q,-D]),sparse(numel(x),2*Etot+2*Bcount)];bv=[zeros(R+2*Bcount,1);-ones(numel(x),1)];
cost=pi*h*[1;areas.*sizes;zeros(2*Etot+2*Bcount,1)];
opts=optimoptions('linprog','Display','none','ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'MaxIterations',20000);
[v,lpval,ef,lpout]=linprog(cost,Amat,bv,Aeq,beq,zeros(numel(cost),1),[],opts);assert(ef>0,'tensor jump LP failed');v=max(v,0);
[fx,fy]=diamond_points(48,192);FD=kmatrix(fx,fy,h,orbits,gn,gw);fv=v(1)*shortQ(fx,fy,h)+FD*v(2:J+1);[dm,idx]=min(fv);assert(dm>0);v=v/dm;
actual_jump=JD*v(1:J+1);jump_parts=v(J+2:J+1+Etot)+v(J+2+Etot:J+1+2*Etot);
out=struct('status','completed_numeric_nonuniform_tensor_jump_PROPOSAL','matlab_version',version,'requested_comp_threads',1,'h',h,'subdiv',subdiv,'grid_id',raw.id,'grid_input_path',grid_input,'axis_edges_exact',{cellstr(tokens)},'axis_intervals',m,'axis_centers_float',center,'axis_widths_float',width,'source_margin',source_margin,'domain_halfwidth',ed(end),'D4_orbit_count',J,'physical_cells',m*m,'source_subrectangles',Bcount,'source_rows',R,'curvature_rows',2*Bcount,'jump_equality_rows',Etot,'LP_variables',numel(v),'LP_nnz',nnz(Amat)+nnz(Aeq),'source_origin_mass',v(1),'density_coefficients',v(2:J+1),'orbit_representatives',rep,'orbit_widths',repwidth,'orbit_cell_areas',areas,'orbit_sizes',sizes,'LP_objective_before_D_rescale',lpval,'LP_exitflag',ef,'LP_output',lpout,'D_coarse_count',numel(x),'D_fine_count',numel(fx),'D_min_before_rescale',dm,'D_minimum_point',[fx(idx),fy(idx)],'cost_from_exact_mass_formula_float',cost'*v,'source_constraint_min_float',-max(sourceA*v),'curvature_constraint_min_float',-max(curvA*v),'jump_equality_max_error_float',max(abs(Aeq*v)),'density_jump_l1_float',sum(abs(actual_jump)),'jump_part_extra_float',sum(jump_parts-abs(actual_jump)),'Qcell_gauss_order',24,'elapsed_seconds',toc(started),'scope','Nonuniform tensor rectangles with center/width D4 symmetry, actual-area mass and row/column-width jumps. Floating proposal only; no general upper without exact source/full-D certificates.');
fid=fopen(fullfile(output_root,'occupation-tensor-jump-lp.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
fprintf('TENSOR-JUMP h=%.6g m=%d s=%d J=%d cells=%d vars=%d cost=%.15g source=%.8g curvature=%.8g equality=%.5g elapsed=%.4g\n',h,m,subdiv,J,m*m,numel(v),cost'*v,out.source_constraint_min_float,out.curvature_constraint_min_float,out.jump_equality_max_error_float,toc(started));
end
function x=fraction_value(token)
p=split(string(token),'/');assert(numel(p)<=2);x=str2double(p(1));if numel(p)==2,d=str2double(p(2));assert(d~=0);x=x/d;end;assert(isfinite(x));
end
function [p,dd]=cauchyjet(x,h)
dmin=min(abs(x));if x(1)<=0&&x(2)>=0,dmin=0;end;dmax=max(abs(x));p=[h/(pi*(h*h+dmax*dmax)),h/(pi*(h*h+dmin*dmin))];
z=[x(1),x(2)];for c=[-h,0,h],if x(1)<=c&&c<=x(2),z(end+1)=c;end,end %#ok<AGROW>
v=2*h*(3*z.^2-h*h)./(pi*(h*h+z.^2).^3);dd=[min(v),max(v)];
end
function [lo,hi]=kprimevec(xl,xr,h)
a=-2*h*xl./(pi*(h*h+xl.^2).^2);b=-2*h*xr./(pi*(h*h+xr.^2).^2);lo=min(a,b);hi=max(a,b);c=h/sqrt(3);mag=9/(8*sqrt(3)*pi*h*h);mask=xl<=c&c<=xr;lo(mask)=-mag;mask=xl<=-c&-c<=xr;hi(mask)=mag;
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
p=atan2(s.*delta,s*s+dist.^2-delta.^2/4)/pi;
end
function D=kmatrix(x,y,h,orbits,gn,gw)
D=zeros(numel(x),numel(orbits));
for j=1:numel(orbits)
 cells=orbits(j).cells;v=zeros(size(x));
 for k=1:size(cells,1)
  dx=abs(x-cells(k,1));dy=abs(y-cells(k,2));
  for q=1:numel(gn),s=h*gn(q);v=v+2*pi*h*gw(q)*boxprob(dx,cells(k,3),s).*boxprob(dy,cells(k,4),s);end
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
