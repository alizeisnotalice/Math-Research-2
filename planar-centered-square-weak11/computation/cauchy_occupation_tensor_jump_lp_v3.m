function cauchy_occupation_tensor_jump_lp_v3(h,grid_input,subdiv,output_root,run_mode,export_sparse)
% R65 V3: canonical orbit-density-difference groups; V1/V2 unchanged.
% Ordinary double LP PROPOSAL only. Root owns launch/resource/exit control.
if nargin<5,run_mode='solve';end;if nargin<6,export_sparse=false;end;assert(ismember(run_mode,{'solve','diagnose','lambda-only'}));assert(isscalar(export_sparse)&&ismember(export_sparse,[false,true]));
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
lambdaDiag=lambda_only_diagnostic(h,subdiv,rep,repwidth);lambdaDiag.grid_id=raw.id;lambdaDiag.axis_edges_exact=cellstr(tokens);lambdaDiag.run_mode=run_mode;
write_json(fullfile(output_root,'tensor-lambda-only-diagnostic.json'),lambdaDiag);
if strcmp(run_mode,'lambda-only'),fprintf('TENSOR-LAMBDA-ONLY h=%.6g minrow=%.12g qmin=%.12g lambda=%.12g numerical_sufficient_feasible=%d\n',h,lambdaDiag.source_row_min_unit_lambda,lambdaDiag.discrete_D_qmin,lambdaDiag.discrete_D_lambda_witness,lambdaDiag.numeric_full_sampled_LP_witness_exists);return;end
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
% Canonical groups use identical density differences, not geometric proximity.
nonzero_physical=find(any(JD(:,2:end),2));pairall=zeros(Etot,2);physical_sign=zeros(Etot,1);
for r=nonzero_physical'
 [~,cols,vals]=find(JD(r,2:end));assert(numel(cols)<=2&&all(abs(vals)==1));
 if numel(cols)==1,pairall(r,:)=[0,cols(1)];physical_sign(r)=vals(1);
 else,pairall(r,:)=[min(cols),max(cols)];physical_sign(r)=full(JD(r,1+max(cols)));end
end
[group_pairs,~,compressed]=unique(pairall(nonzero_physical,:),'rows');G=size(group_pairs,1);physical_group=zeros(Etot,1);physical_group(nonzero_physical)=compressed;
GD=spalloc(G,J+1,2*G);
for g=1:G
 lo=group_pairs(g,1);hi=group_pairs(g,2);assert(hi>lo&&hi<=J);GD(g,1+hi)=1;if lo>0,GD(g,1+lo)=-1;end
end
SG=sparse(nonzero_physical,physical_group(nonzero_physical),physical_sign(nonzero_physical),Etot,G);
assert(nnz(JD-SG*GD)==0,'canonical grouping changed the exact integer density difference map');
GX=SG(1:Ex,:);GY=SG(Ex+1:end,:);XP=max(GX,0);XN=max(-GX,0);YP=max(GY,0);YN=max(-GY,0);
group_counts=full(sum(abs(SG),1))';zero_interface_count=Etot-numel(nonzero_physical);
Bcount=J*subdiv^2;R=4*Bcount;V=zeros(R,J+1);target=zeros(R,1);ox=zeros(Bcount,1);oy=ox;
PxL=zeros(Bcount,G);PxU=PxL;PyL=zeros(Bcount,G);PyU=PyL;errx=zeros(R,1);erry=errx;bi=0;
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
  [pl,pu]=mulvec(kl,ku,tl,tu);PxL(bi,:)=pl'*XP-pu'*XN;PxU(bi,:)=pu'*XP-pl'*XN;
  [kl,ku]=kprimevec(B(2,1)-ye,B(2,2)-ye,h);[tl,tu]=boxrangevec(B(1,1)-yx,B(1,2)-yx,yxw,h);
  [pl,pu]=mulvec(kl,ku,tl,tu);PyL(bi,:)=pl'*YP-pu'*YN;PyU(bi,:)=pu'*YP-pl'*YN;
 end,end
end
source_margin=1e-7;V(:,1)=V(:,1)-source_margin*max(1,abs(V(:,1)));
baseA=-V;baseA(sub2ind(size(baseA),(1:R)',target+1))=baseA(sub2ind(size(baseA),(1:R)',target+1))+1;
% Anisotropic interpolation errors use each target subrectangle's own lengths.
EX=sparse((1:R)',repelem((1:Bcount)',4),errx,R,Bcount);EY=sparse((1:R)',repelem((1:Bcount)',4),erry,R,Bcount);
sourceA=[sparse(baseA),sparse(R,2*G),EX,EY];
Cx=sparse((1:Bcount)',ones(Bcount,1),ox,Bcount,J+1);Cy=sparse((1:Bcount)',ones(Bcount,1),oy,Bcount,J+1);
ZZ=sparse(Bcount,Bcount);II=speye(Bcount);
curvA=[Cx,sparse(PxU),sparse(-PxL),-II,ZZ;Cy,sparse(PyU),sparse(-PyL),ZZ,-II];
Aeq=[GD,-speye(G),speye(G),sparse(G,2*Bcount)];beq=zeros(G,1);
[gn,gw]=gauss(24);[x,y]=diamond_points(24,96);D=kmatrix(x,y,h,orbits,gn,gw);q=shortQ(x,y,h);
Amat=[sourceA;curvA;sparse([-q,-D]),sparse(numel(x),2*G+2*Bcount)];bv=[zeros(R+2*Bcount,1);-ones(numel(x),1)];
cost=pi*h*[1;areas.*sizes;zeros(2*G+2*Bcount,1)];
opts=optimoptions('linprog','Display','none','ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'MaxIterations',20000);
matrixDiag=struct('physical_interfaces',Etot,'nonzero_physical_interfaces',numel(nonzero_physical),'zero_interface_count',zero_interface_count,'canonical_jump_groups',G,'group_multiplicity_range',[min(group_counts),max(group_counts)],'source_density_orbits',J,'status','BUILT_BEFORE_SOLVER','run_mode',run_mode,'h',h,'grid_id',raw.id,'subdiv',subdiv,'lambda_only',lambdaDiag,'inequality_matrix',matrix_summary(Amat),'equality_matrix',matrix_summary(Aeq),'source_matrix',matrix_summary(sourceA),'curvature_matrix',matrix_summary(curvA),'cost',vector_summary(cost),'inequality_rhs',vector_summary(bv),'q_D',vector_summary(q),'D_matrix',matrix_summary(D),'LP_MaxIterations',20000,'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'build_elapsed_seconds',toc(started),'scope','Floating formulation/scale diagnostic, not infeasibility proof or wholeD certificate.');
write_json(fullfile(output_root,'tensor-matrix-diagnostic.json'),matrixDiag);
if export_sparse
 lower_bounds=zeros(numel(cost),1);save(fullfile(output_root,'tensor-lp-matrices.mat'),'Amat','bv','Aeq','beq','cost','lower_bounds','h','subdiv','ed','center','width','rep','repwidth','sizes','areas','JD','GD','SG','group_pairs','physical_group','physical_sign','J','G','Bcount','R','-v7');
end
if strcmp(run_mode,'diagnose'),fprintf('TENSOR-DIAGNOSE matrices built h=%.6g minlambdaonlyrow=%.10g elapsed=%.4g\n',h,lambdaDiag.source_row_min_unit_lambda,toc(started));return;end
v=[];lpval=[];ef=[];lpout=struct();solverException=struct();solverTimer=tic;
try
 [v,lpval,ef,lpout]=linprog(cost,Amat,bv,Aeq,beq,zeros(numel(cost),1),[],opts);
catch err
 solverException=struct('identifier',err.identifier,'message',err.message,'stack',err.stack,'report',getReport(err,'extended','hyperlinks','off'));
end
solverStatus=struct('status','LP_RETURNED','h',h,'grid_id',raw.id,'subdiv',subdiv,'run_mode',run_mode,'exitflag',ef,'objective',lpval,'LP_output',lpout,'exception',solverException,'returned_vector',v,'solver_seconds',toc(solverTimer),'elapsed_seconds',toc(started),'matrix_diagnostic_file','tensor-matrix-diagnostic.json','lambda_only',lambdaDiag,'scope','Solver status only; failure does not prove mathematical kernel or even sufficient-row infeasibility.');
if ~isempty(v)&&numel(v)==numel(cost)&&all(isfinite(v)),solverStatus.inequality_max_violation=max(Amat*v-bv);solverStatus.equality_max_abs_residual=max(abs(Aeq*v-beq));solverStatus.lower_bound_violation=max(-v);end
if isempty(ef)||ef<=0,solverStatus.status='LP_FAILED_DIAGNOSTIC_SAVED';end
write_json(fullfile(output_root,'tensor-lp-solver-status.json'),solverStatus);
if isempty(ef)||ef<=0,error('R65:LPFailedAfterDiagnostic','LP did not succeed; full exitflag/output in tensor-lp-solver-status.json');end
v=max(v,0);
[fx,fy]=diamond_points(48,192);FD=kmatrix(fx,fy,h,orbits,gn,gw);fv=v(1)*shortQ(fx,fy,h)+FD*v(2:J+1);[dm,idx]=min(fv);assert(dm>0);v=v/dm;
actual_jump=GD*v(1:J+1);jump_parts=v(J+2:J+1+G)+v(J+2+G:J+1+2*G);
out=struct('status','completed_numeric_nonuniform_tensor_jump_V3_GROUPED_PROPOSAL','matlab_version',version,'requested_comp_threads',1,'h',h,'subdiv',subdiv,'grid_id',raw.id,'grid_input_path',grid_input,'axis_edges_exact',{cellstr(tokens)},'axis_intervals',m,'axis_centers_float',center,'axis_widths_float',width,'source_margin',source_margin,'domain_halfwidth',ed(end),'D4_orbit_count',J,'physical_cells',m*m,'source_subrectangles',Bcount,'source_rows',R,'curvature_rows',2*Bcount,'physical_interfaces',Etot,'zero_interface_count',zero_interface_count,'canonical_jump_groups',G,'group_multiplicity_range',[min(group_counts),max(group_counts)],'jump_equality_rows',G,'LP_variables',numel(v),'LP_nnz',nnz(Amat)+nnz(Aeq),'source_origin_mass',v(1),'density_coefficients',v(2:J+1),'orbit_representatives',rep,'orbit_widths',repwidth,'orbit_cell_areas',areas,'orbit_sizes',sizes,'LP_objective_before_D_rescale',lpval,'LP_exitflag',ef,'LP_output',lpout,'D_coarse_count',numel(x),'D_fine_count',numel(fx),'D_min_before_rescale',dm,'D_minimum_point',[fx(idx),fy(idx)],'cost_from_exact_mass_formula_float',cost'*v,'source_constraint_min_float',-max(sourceA*v),'curvature_constraint_min_float',-max(curvA*v),'jump_equality_max_error_float',max(abs(Aeq*v)),'group_density_jump_l1_float',sum(abs(actual_jump)),'physical_density_jump_l1_float',sum(abs(JD*v(1:J+1))),'jump_part_extra_float',sum(jump_parts-abs(actual_jump)),'Qcell_gauss_order',24,'elapsed_seconds',toc(started),'scope','Canonical density-difference shared groups with signed physical-interface interval sums, nonuniform tensor actual-area mass. Exact function identity; floating sufficient LP may differ numerically. No general upper without source/full-D certificates.');
fid=fopen(fullfile(output_root,'occupation-tensor-jump-lp.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
fprintf('TENSOR-JUMP h=%.6g m=%d s=%d J=%d cells=%d vars=%d cost=%.15g source=%.8g curvature=%.8g equality=%.5g elapsed=%.4g\n',h,m,subdiv,J,m*m,numel(v),cost'*v,out.source_constraint_min_float,out.curvature_constraint_min_float,out.jump_equality_max_error_float,toc(started));
end
function out=lambda_only_diagnostic(h,subdiv,rep,repwidth)
% Deterministic sufficient-row witness using f=0; NO LP/no rectangle basis.
J=size(rep,1);Bcount=J*subdiv^2;minrow=inf;worst=struct();bi=0;Mxx=zeros(Bcount,1);Myy=Mxx;source_margin=1e-7;
for i=1:J
 wx=repwidth(i,1);wy=repwidth(i,2);
 for ix=0:subdiv-1,for iy=0:subdiv-1
  bi=bi+1;B=[rep(i,1)-wx/2+wx*[ix,ix+1]/subdiv;rep(i,2)-wy/2+wy*[iy,iy+1]/subdiv];
  vertices=[B(1,1),B(2,1);B(1,1),B(2,2);B(1,2),B(2,1);B(1,2),B(2,2)];
  pv=ph(vertices(:,1),vertices(:,2),h);[px,pxx]=cauchyjet(B(1,:),h);[py,pyy]=cauchyjet(B(2,:),h);ux=imul(pxx,py);uy=imul(px,pyy);
  Mxx(bi)=max(0,ux(2));Myy(bi)=max(0,uy(2));ex=(wx/subdiv)^2*Mxx(bi)/8;ey=(wy/subdiv)^2*Myy(bi)/8;
  rows=pv-source_margin*max(1,abs(pv))-ex-ey;[v,k]=min(rows);
  if v<minrow,minrow=v;worst=struct('target_orbit_1based',i,'subrectangle_1based',bi,'target_center',rep(i,:),'target_width',repwidth(i,:),'box',B,'vertex',vertices(k,:),'origin_p_vertex',pv(k),'origin_curvature_upper',[ux(2),uy(2)],'eta_unit_lambda',[Mxx(bi),Myy(bi)],'interpolation_errors',[ex,ey],'source_row_coefficient',v);end
 end,end
end
[x,y]=diamond_points(24,96);q=shortQ(x,y,h);[qmin,qi]=min(q);qmax=max(q);lambda=[];numeric=false;
if isfinite(qmin)&&qmin>0,lambda=(1+1e-6)/qmin;numeric=minrow>=0;end
out=struct('status','FLOAT_LAMBDA_ONLY_SUFFICIENT_ROWS','h',h,'subdiv',subdiv,'D4_orbits',J,'target_subrectangles',Bcount,'source_rows',4*Bcount,'source_margin',source_margin,'source_row_min_unit_lambda',minrow,'worst_row',worst,'discrete_D_qmin',qmin,'discrete_D_qmax',qmax,'discrete_D_qmin_point',[x(qi),y(qi)],'discrete_D_lambda_witness',lambda,'numeric_full_sampled_LP_witness_exists',numeric,'unit_lambda_eta_x',Mxx,'unit_lambda_eta_y',Myy,'scope','f=0 gives a mathematically legal source lambda*p_h independently of row diagnostic. Floating sufficient-row/D sampling feasibility only, not fullD or infeasibility proof.');
if ~isempty(lambda),out.lambda_only_sampled_cost=pi*h*lambda;out.source_row_min_scaled=lambda*minrow;out.D_sampled_margin=1e-6;end
end
function out=matrix_summary(A)
v=nonzeros(A);r=full(max(abs(A),[],2));c=full(max(abs(A),[],1));
out=struct('rows',size(A,1),'columns',size(A,2),'nnz',nnz(A),'zero_rows',sum(r==0),'zero_columns',sum(c==0),'row_inf_range',[min(r),max(r)],'column_inf_range',[min(c),max(c)]);
if ~isempty(v),out.nonzero_abs_range=[min(abs(v)),max(abs(v))];out.signed_range=[min(v),max(v)];out.positive_nnz=sum(v>0);out.negative_nnz=sum(v<0);end
end
function out=vector_summary(v)
v=v(:);out=struct('count',numel(v),'finite_count',sum(isfinite(v)),'zero_count',sum(v==0));
if ~isempty(v),out.range=[min(v),max(v)];z=v(v~=0);if ~isempty(z),out.nonzero_abs_range=[min(abs(z)),max(abs(z))];end,end
end
function write_json(path,out)
fid=fopen(path,'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
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
