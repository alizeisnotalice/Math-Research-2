function cauchy_occupation_cell_lp(h,delta,N,mode,output_root)
% R60 positive occupation kernel proposal; all arithmetic is nondirected.
assert(ismember(mode,{'midpoint','cellinf'}));assert(h>0&&delta>0&&N>=1);
assert(startsWith(output_root,'/Volumes/data/'));if ~isfolder(output_root),mkdir(output_root);end
started=tic;orbits=struct('center',{},'points',{});rep=[];J=0;
for i=0:N,for j=0:i
 z=delta*[i,j];pts=[];for sx=[-1,1],for sy=[-1,1],pts=[pts;sx*z(1),sy*z(2);sx*z(2),sy*z(1)];end,end %#ok<AGROW>
 J=J+1;orbits(J).center=z;orbits(J).points=unique(pts,'rows');rep=[rep;z]; %#ok<AGROW>
end,end
sizes=arrayfun(@(o)size(o.points,1),orbits)';L=zeros(J,J);b=ph(rep(:,1),rep(:,2),h);
if strcmp(mode,'cellinf'),b=ph(rep(:,1)+delta/2,rep(:,2)+delta/2,h);end
for j=1:J
 pts=orbits(j).points;
 for k=1:size(pts,1)
  dx=abs(rep(:,1)-pts(k,1));dy=abs(rep(:,2)-pts(k,2));
  if strcmp(mode,'cellinf'),dx=dx+delta/2;dy=dy+delta/2;end
  L(:,j)=L(:,j)+boxprob(dx,delta,h).*boxprob(dy,delta,h);
 end
end
[gn,gw]=gauss(24);[x,y]=diamond_points(32,128);D=kmatrix(x,y,h,delta,orbits,gn,gw);q=shortQ(x,y,h);
Amat=[-b,speye(J)-sparse(L);-q,-D];bv=[zeros(J,1);-ones(numel(x),1)];cost=pi*h*[1;delta^2*sizes];opts=optimoptions('linprog','Display','none','ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
[v,lpval,ef,lpout]=linprog(cost,Amat,bv,[],[],zeros(J+1,1),[],opts);assert(ef>0,'occupation LP failed');v=max(v,0);
[fx,fy]=diamond_points(64,256);FD=kmatrix(fx,fy,h,delta,orbits,gn,gw);fv=v(1)*shortQ(fx,fy,h)+FD*v(2:end);[dm,idx]=min(fv);assert(dm>0);v=v/dm;
% Actual source residual sampled at every representative cell's center and corners.
% In cellinf mode entrywise lower constraints are sufficient analytically, but still float here.
shift=[0,0;-.5,-.5;-.5,.5;.5,-.5;.5,.5]*delta;
resmin=inf;respoint=[];rescell=0;
for qi=1:size(shift,1)
 xx=rep(:,1)+shift(qi,1);yy=rep(:,2)+shift(qi,2);rr=v(1)*ph(xx,yy,h)-v(2:end);
 for j=1:J
  if v(j+1)==0,continue;end
  pts=orbits(j).points;
  for k=1:size(pts,1),rr=rr+v(j+1)*boxprob(abs(xx-pts(k,1)),delta,h).*boxprob(abs(yy-pts(k,2)),delta,h);end
 end
 [rmin,ri]=min(rr);if rmin<resmin,resmin=rmin;respoint=[xx(ri),yy(ri)];rescell=ri;end
end
lowerres=v(1)*b+L*v(2:end)-v(2:end);
out=struct('status','completed_numeric_positive_occupation_LP','matlab_version',version,'threads',1,'h',h,'delta',delta,'N',N,'mode',mode,'domain_halfwidth',(N+.5)*delta,'D4_orbit_count',J,'physical_cells',sum(sizes),'source_origin_mass',v(1),'density_coefficients',v(2:end),'orbit_representatives',rep,'orbit_sizes',sizes,'cell_transition_row_sum_range',[min(sum(L,2)),max(sum(L,2))],'LP_objective_before_D_rescale',lpval,'LP_exitflag',ef,'LP_output',lpout,'D_coarse_count',numel(x),'D_fine_count',numel(fx),'D_min_before_rescale',dm,'D_minimum_point',[fx(idx),fy(idx)],'cost_from_exact_mass_formula_float',cost'*v,'source_residual_constraint_min_float',min(lowerres),'source_residual_sample_min',resmin,'source_residual_sample_point',respoint,'source_residual_cell_index',rescell,'Qcell_gauss_order',24,'elapsed_seconds',toc(started),'scope','Positive K=Qh*(A delta0+piecewiseconstant f), exact analytical halfmass=pi h(A+integral f). Midpoint mode does NOT ensure nonnegative source. Cellinf inequalities are analytically sufficient but float LP/entries and D samples are NOT certified. No general upper until directed source positivity/fullD certificate.');
fid=fopen(fullfile(output_root,'occupation-cell-lp.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
fprintf('OCC h=%.6g delta=%.6g N=%d mode=%s J=%d LP=%.15g cost=%.15g Dmin=%.10g residual=%.8g elapsed=%.4g\n',h,delta,N,mode,J,lpval,cost'*v,dm,resmin,toc(started));
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
